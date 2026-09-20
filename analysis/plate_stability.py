import json, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from analyze_glds7 import read_matrix, read_annotation, clean, bh

ROOT=Path('/workspace/scratch/1b7c2e5453f2'); SITE=ROOT/'orbital-cell-twin'

def main():
    expr=pd.read_csv(ROOT/'GSE56659_RMA_expression.csv.gz',index_col=0)
    _,meta=read_matrix(); ann=read_annotation().set_index('ID')
    loci=ann['Platform_ORF'].reindex(expr.index).map(clean).str.upper();valid=loci.str.match(r'^AT[1-5CM]G\d{5}$',na=False)
    pick=pd.DataFrame({'probe':expr.index[valid],'locus':loci[valid].values,'v':expr.loc[valid].var(axis=1).values}).sort_values('v',ascending=False).drop_duplicates('locus')
    G=expr.loc[pick.probe];G.index=pick.locus.values
    mapping={};names={}
    for probe,row in ann.iterrows():
        locus=clean(row.get('Platform_ORF','')).upper()
        if not re.match(r'^AT[1-5CM]G\d{5}$',locus):continue
        ids=[x.strip() for x in str(row.get('GO:Process ID','')).split('///') if x.strip() and x!='nan']
        labs=[x.strip() for x in str(row.get('GO:Process','')).split('///') if x.strip() and x!='nan']
        for go,lab in zip(ids,labs):mapping.setdefault(go,set()).add(locus);names[go]=lab
    focus=pd.read_csv(SITE/'rma-defense-confirmation.csv');root=meta[meta.tissue=='root'].copy();root['plate']=root.title.str.extract(r'plate(\d)',expand=False)
    rows=[]
    for _,x in focus.iterrows():
        genes=sorted(mapping.get(x.go_id,set())&set(G.index));scores=G.loc[genes,root.accession].mean(axis=0)
        effects=[]
        for plate in ['2','4','6']:
            m=root[root.plate==plate];f=scores[m.loc[m.condition=='flight','accession']].mean();g=scores[m.loc[m.condition=='ground','accession']].mean();effects.append(float(f-g))
        rows.append({'go_id':x.go_id,'term':names.get(x.go_id,x.term),'genes_in_score':len(genes),
                     'plate2_effect':effects[0],'plate4_effect':effects[1],'plate6_effect':effects[2],
                     'positive_plates':sum(e>0 for e in effects),'stable_all_plates':all(e>0 for e in effects),
                     'mean_effect':float(np.mean(effects))})
    out=pd.DataFrame(rows).sort_values(['stable_all_plates','mean_effect'],ascending=[False,False]);out.to_csv(SITE/'plate-stability.csv',index=False)
    payload={'method':'Mean RMA expression of all genes annotated to each selected GO term; flight-minus-ground effect calculated separately within plates 2, 4, and 6.',
             'stable_rule':'Positive flight effect in all three plates','stable_terms':int(out.stable_all_plates.sum()),
             'tested_terms':len(out),'terms':out.to_dict('records'),
             'caution':'Post-selection sensitivity analysis; plates are experimental blocks, not independent space missions.'}
    (SITE/'plate-stability.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(payload,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
