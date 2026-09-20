import json
import re
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

ROOT = Path('/workspace/scratch/1b7c2e5453f2')
sys.path.insert(0, str(ROOT/'pyaffy-src'))
import pyaffy

from analyze_glds7 import read_matrix, read_annotation, clean, bh

RAW = ROOT/'GSE56659_RAW'
CDF = ROOT/'ATH1-121501.CDF'
SITE = ROOT/'orbital-cell-twin'


def prior(s2, df):
    z=np.log(np.maximum(s2,np.finfo(float).tiny))-special.digamma(df/2)+np.log(df/2)
    target=max(float(np.var(z,ddof=1)-special.polygamma(1,df/2)),1e-8)
    try: df0=float(optimize.brentq(lambda x:special.polygamma(1,x/2)-target,1e-3,1e7))
    except ValueError: df0=1e7
    s02=float(np.exp(np.mean(z)-special.digamma(df0/2)+np.log(df0/2)))
    return df0,s02


def design_test(Y, meta):
    # condition plus plate blocking; intercept is explicit
    plate=meta.title.str.extract(r'plate(\d)',expand=False)
    X=pd.DataFrame({'intercept':np.ones(len(meta)),'flight':(meta.condition=='flight').astype(float).to_numpy(),
                    'plate4':(plate=='4').astype(float).to_numpy(),'plate6':(plate=='6').astype(float).to_numpy()},index=meta.accession)
    X=X.loc[:,X.var().gt(0)|X.columns.to_series().eq('intercept')]
    A=X.to_numpy(float); inv=np.linalg.pinv(A.T@A); beta=inv@A.T@Y.T
    fitted=A@beta; resid=Y.T-fitted; df=A.shape[0]-np.linalg.matrix_rank(A)
    s2=(resid**2).sum(axis=0)/df; df0,s02=prior(s2,df); post=(df*s2+df0*s02)/(df+df0)
    j=list(X.columns).index('flight'); se=np.sqrt(post*inv[j,j]); t=beta[j]/se
    p=2*stats.t.sf(np.abs(t),df+df0)
    return beta[j],t,p,df,df0,s02


def go_members(ann):
    out={}
    for probe,row in ann.iterrows():
        locus=clean(row.get('Platform_ORF','')).upper()
        if not re.match(r'^AT[1-5CM]G\d{5}$',locus): continue
        ids=[x.strip() for x in str(row.get('GO:Process ID','')).split('///') if x.strip() and x.strip()!='nan']
        names=[x.strip() for x in str(row.get('GO:Process','')).split('///') if x.strip() and x.strip()!='nan']
        for go,name in zip(ids,names): out.setdefault(go,{'term':name,'genes':set()})['genes'].add(locus)
    return out


def main():
    _,meta=read_matrix(); files={p.name.split('_')[0]:p for p in RAW.glob('*.CEL.gz')}
    sample_files=OrderedDict((a,str(files[a])) for a in meta.accession)
    genes,samples,X=pyaffy.rma(str(CDF),sample_files)
    probes=[]
    for g in genes:
        value=g.decode() if isinstance(g,bytes) else str(g)
        if value.startswith("b'") and value.endswith("'"):
            value=value[2:-1]
        probes.append(value)
    expr=pd.DataFrame(X,index=probes,columns=samples)
    expr.to_csv(ROOT/'GSE56659_RMA_expression.csv.gz',compression='gzip')
    ann=read_annotation().set_index('ID'); loci=ann['Platform_ORF'].reindex(expr.index).map(clean).str.upper()
    valid=loci.str.match(r'^AT[1-5CM]G\d{5}$',na=False)
    # representative probe selected by overall variance, independent of condition
    pick=pd.DataFrame({'probe':expr.index[valid],'locus':loci[valid].values,'v':expr.loc[valid].var(axis=1).values}).sort_values('v',ascending=False).drop_duplicates('locus')
    G=expr.loc[pick.probe]; G.index=pick.locus.values
    rows=[]; priors={}
    for tissue in ['whole plant','shoot','hypocotyl','root']:
        m=meta[meta.tissue==tissue].copy(); Y=G[m.accession].to_numpy()
        fc,t,p,df,df0,s02=design_test(Y,m); q=bh(p)
        part=pd.DataFrame({'tissue':tissue,'locus':G.index,'log2fc':fc,'moderated_t':t,'p':p,'fdr':q})
        rows.append(part); priors[tissue]={'residual_df':int(df),'prior_df':df0,'prior_variance':s02}
    res=pd.concat(rows,ignore_index=True);res.to_csv(SITE/'rma-moderated-results.csv',index=False)

    # Confirm the previously identified root defense terms using the independent RMA preprocessing.
    go=go_members(ann); root=res[res.tissue=='root']; universe=set(root.locus)
    focus=pd.read_csv(SITE/'go-enrichment.csv');focus=focus[(focus.tissue=='root')&(focus.direction=='up')&(focus.fdr<.05)].sort_values('fdr').head(12)
    selected=set(root.loc[(root.p<=.01)&(root.log2fc>0),'locus']); out=[]
    for _,x in focus.iterrows():
        members=go.get(x.go_id,{'genes':set()})['genes']&universe; overlap=selected&members
        pv=stats.hypergeom.sf(len(overlap)-1,len(universe),len(members),len(selected)) if overlap else 1.0
        out.append({'go_id':x.go_id,'term':go.get(x.go_id,{'term':x.term})['term'],'selected_genes':len(selected),
                    'term_size':len(members),'overlap':len(overlap),'p':pv,'genes':','.join(sorted(overlap))})
    q=bh([x['p'] for x in out])
    for x,z in zip(out,q):x['fdr']=float(z);x['confirmed']=bool(z<.05 and x['overlap']>=2)
    pd.DataFrame(out).sort_values('fdr').to_csv(SITE/'rma-defense-confirmation.csv',index=False)

    # Pipeline sensitivity: rank agreement of RMA and original MAS5 effects.
    old=pd.read_csv(SITE/'moderated-gene-results.csv'); comp=[]
    for tissue in ['whole plant','shoot','hypocotyl','root']:
        a=res[res.tissue==tissue].set_index('locus');b=old[old.tissue==tissue].set_index('locus');idx=a.index.intersection(b.index)
        rho,pv=stats.spearmanr(a.loc[idx,'log2fc'],b.loc[idx,'log2fc'])
        top=set(a.nsmallest(300,'p').index);top_old=set(b.nsmallest(300,'p').index);inter=top&top_old
        same=sum(np.sign(a.loc[g,'log2fc'])==np.sign(b.loc[g,'log2fc']) for g in inter)
        comp.append({'tissue':tissue,'shared_loci':len(idx),'effect_spearman':rho,'spearman_p':pv,
                     'top300_overlap':len(inter),'top300_same_direction':same/len(inter) if inter else None})
    pd.DataFrame(comp).to_csv(SITE/'pipeline-sensitivity.csv',index=False)
    payload={'raw_data':{'source':'GEO GSE56659_RAW.tar','cel_files':len(files),'cdf':'ATH1-121501.CDF',
                         'preprocessing':'RMA background correction, quantile normalization, median polish'},
             'model':'Tissue-specific flight effect blocked by plate (2, 4, 6), empirical-Bayes variance moderation',
             'genes_tested':int(G.shape[0]),'priors':priors,
             'tissue_summary':[{'tissue':t,'fdr_005':int((g.fdr<.05).sum()),'p001':int((g.p<=.01).sum())} for t,g in res.groupby('tissue',sort=False)],
             'defense_confirmation':{'tested_terms':len(out),'confirmed_terms':sum(x['confirmed'] for x in out),'terms':out},
             'pipeline_sensitivity':comp,
             'caution':'RMA was reproduced in Python from raw CEL files; confirmatory publication should independently rerun the frozen workflow in Bioconductor and archive checksums.'}
    (SITE/'raw-rma-confirmation.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(payload,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
