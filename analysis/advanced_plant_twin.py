import json
import numpy as np
import pandas as pd
from scipy import optimize, special, stats
from sklearn.metrics import roc_auc_score

from analyze_glds7 import read_matrix, read_annotation, clean, bh
from cross_validate_nasa import external_sets
from config import RESULTS_DIR


SITE = RESULTS_DIR
JSON_OUT = SITE / 'advanced-validation.json'
GENE_OUT = SITE / 'moderated-gene-results.csv'
GO_OUT = SITE / 'go-enrichment.csv'
TRANSFER_OUT = SITE / 'transfer-validation.csv'


def estimate_prior_variance(s2, df_resid):
    """Estimate scaled inverse-chi-square prior from residual variances."""
    s2 = np.asarray(s2, float)
    s2 = np.maximum(s2[np.isfinite(s2)], np.finfo(float).tiny)
    z = np.log(s2) - special.digamma(df_resid / 2) + np.log(df_resid / 2)
    target = max(float(np.var(z, ddof=1) - special.polygamma(1, df_resid / 2)), 1e-8)
    f = lambda x: special.polygamma(1, x / 2) - target
    try:
        df0 = float(optimize.brentq(f, 1e-3, 1e7))
    except ValueError:
        df0 = 1e7
    log_s02 = float(np.mean(z) - special.digamma(df0 / 2) + np.log(df0 / 2))
    return df0, float(np.exp(log_s02))


def split_terms(value):
    if pd.isna(value) or not str(value).strip() or str(value).strip() == '---':
        return []
    return [x.strip() for x in str(value).split('///') if x.strip() and x.strip() != '---']


def prepare_gene_data():
    expr, meta = read_matrix()
    ann = read_annotation().set_index('ID')
    log_expr = np.log2(expr + 1)
    loci = ann['Platform_ORF'].reindex(expr.index).map(clean).str.upper()
    valid = loci.str.match(r'^AT[1-5CM]G\d{5}$', na=False)
    # One probe per locus, selected by overall variance without using group labels.
    probe_map = pd.DataFrame({'probe': expr.index[valid], 'locus': loci[valid].values,
                              'variance': log_expr.loc[valid].var(axis=1).values})
    chosen = probe_map.sort_values('variance', ascending=False).drop_duplicates('locus')
    genes = log_expr.loc[chosen.probe].copy()
    genes.index = chosen.locus.values
    probe_for_gene = dict(zip(chosen.locus, chosen.probe))
    return genes, meta, ann, probe_for_gene


def moderated_results(genes, meta, probe_for_gene):
    frames, priors = [], {}
    for tissue in ['whole plant', 'shoot', 'hypocotyl', 'root']:
        fids = meta.loc[(meta.tissue == tissue) & (meta.condition == 'flight'), 'accession'].tolist()
        gids = meta.loc[(meta.tissue == tissue) & (meta.condition == 'ground'), 'accession'].tolist()
        f, g = genes[fids].to_numpy(), genes[gids].to_numpy()
        lfc = f.mean(axis=1) - g.mean(axis=1)
        df = len(fids) + len(gids) - 2
        s2 = ((len(fids)-1)*f.var(axis=1, ddof=1) + (len(gids)-1)*g.var(axis=1, ddof=1)) / df
        df0, s02 = estimate_prior_variance(s2, df)
        post = (df*s2 + df0*s02) / (df + df0)
        se = np.sqrt(post * (1/len(fids) + 1/len(gids)))
        t = np.divide(lfc, se, out=np.zeros_like(lfc), where=se > 0)
        p = 2 * stats.t.sf(np.abs(t), df + df0)
        frame = pd.DataFrame({'tissue': tissue, 'locus': genes.index,
                              'probe': [probe_for_gene[x] for x in genes.index],
                              'log2fc': lfc, 'moderated_t': t, 'p': p})
        frame['fdr'] = bh(frame.p)
        frames.append(frame)
        priors[tissue] = {'residual_df': df, 'prior_df': df0, 'prior_variance': s02,
                          'n_flight': len(fids), 'n_ground': len(gids)}
    return pd.concat(frames, ignore_index=True), priors


def go_enrichment(gene_res, ann, probe_for_gene):
    gene_to_terms, term_name = {}, {}
    for locus, probe in probe_for_gene.items():
        if probe not in ann.index:
            continue
        names = split_terms(ann.at[probe, 'GO:Process'])
        ids = split_terms(ann.at[probe, 'GO:Process ID'])
        pairs = list(zip(ids, names)) if len(ids) == len(names) else [(x, x) for x in ids]
        gene_to_terms[locus] = {i for i, _ in pairs}
        term_name.update({i: n for i, n in pairs})
    universe = set(gene_to_terms)
    term_to_genes = {}
    for gene, terms in gene_to_terms.items():
        for term in terms:
            term_to_genes.setdefault(term, set()).add(gene)
    term_to_genes = {k:v for k,v in term_to_genes.items() if 5 <= len(v) <= 500}
    rows = []
    for tissue, sub in gene_res.groupby('tissue'):
        for direction in ['up', 'down']:
            selected = set(sub.loc[(sub.p <= .01) & ((sub.log2fc > 0) if direction == 'up' else (sub.log2fc < 0)), 'locus']) & universe
            if not selected:
                continue
            local = []
            for term, members in term_to_genes.items():
                overlap = selected & members
                if not overlap:
                    continue
                p = stats.hypergeom.sf(len(overlap)-1, len(universe), len(members), len(selected))
                local.append({'tissue':tissue, 'direction':direction, 'go_id':term,
                              'term':term_name.get(term, term), 'selected_genes':len(selected),
                              'term_size':len(members), 'overlap':len(overlap), 'p':p,
                              'genes':','.join(sorted(overlap)[:25])})
            if local:
                q = bh([x['p'] for x in local])
                for x, fdr in zip(local, q): x['fdr'] = fdr
                rows.extend(local)
    return pd.DataFrame(rows).sort_values(['fdr','p']) if rows else pd.DataFrame()


def transfer_validation(gene_res):
    ext = external_sets()
    rows=[]
    for tissue, a in gene_res.groupby('tissue'):
        a=a.set_index('locus')
        for study,b in ext.items():
            common=a.index.intersection(b.index)
            av=a.loc[common,'log2fc'].astype(float)
            bv=b.loc[common,'log2fc'].astype(float)
            rho,rho_p=stats.spearmanr(av,bv)
            sig=b.loc[common][b.loc[common,'p']<=.01]
            labels=(sig.log2fc>0).astype(int)
            scores=a.loc[sig.index,'log2fc']
            if len(labels)>=10 and labels.nunique()==2:
                auc=roc_auc_score(labels,scores)
                u_p=stats.mannwhitneyu(scores[labels==1],scores[labels==0],alternative='two-sided').pvalue
                acc=float((np.sign(scores)==np.sign(sig.log2fc)).mean())
            else:
                auc=acc=np.nan; u_p=1.0
            rows.append({'tissue':tissue,'study':study,'shared_genes':len(common),
                         'spearman_rho':rho,'spearman_p':rho_p,
                         'external_p001_genes':len(sig),'direction_auc':auc,
                         'direction_accuracy':acc,'auc_p':u_p})
    frame=pd.DataFrame(rows)
    frame['spearman_fdr']=bh(frame.spearman_p.fillna(1))
    frame['auc_fdr']=bh(frame.auc_p.fillna(1))
    frame['transfer_status']=np.where((frame.auc_fdr<.05)&(frame.direction_auc>=.6),'predictive',
                              np.where((frame.spearman_fdr<.05)&(frame.spearman_rho>0),'rank-consistent','not supported'))
    return frame.sort_values(['transfer_status','auc_fdr','spearman_fdr'])


def main():
    genes, meta, ann, probe_for_gene = prepare_gene_data()
    gene_res, priors = moderated_results(genes, meta, probe_for_gene)
    go = go_enrichment(gene_res, ann, probe_for_gene)
    transfer = transfer_validation(gene_res)
    gene_res.to_csv(GENE_OUT,index=False)
    go.to_csv(GO_OUT,index=False)
    transfer.to_csv(TRANSFER_OUT,index=False)
    tissue_summary=[]
    for tissue,sub in gene_res.groupby('tissue',sort=False):
        tissue_summary.append({'tissue':tissue,'genes_tested':len(sub),
          'fdr_005':int((sub.fdr<.05).sum()),'p001':int((sub.p<=.01).sum()),
          'abs_fc1_p005':int(((sub.p<=.05)&(sub.log2fc.abs()>=1)).sum())})
    top_go=[] if go.empty else go.head(15).replace({np.nan:None}).to_dict('records')
    supported = transfer[transfer.transfer_status != 'not supported'].copy()
    supported['priority'] = supported.transfer_status.map({'predictive':0, 'rank-consistent':1})
    best_transfer=supported.sort_values(['priority','auc_fdr','spearman_fdr']).drop(columns='priority').head(12).replace({np.nan:None}).to_dict('records')
    payload={'method':{
      'expression':'log2(MAS5 intensity + 1)','unit':'one highest-variance probe per Arabidopsis locus',
      'differential_model':'two-group linear contrast by tissue with empirical-Bayes moderated residual variance',
      'multiple_testing':'Benjamini-Hochberg within tissue; GO within tissue and direction; transfer across 20 comparisons',
      'go':'one-sided hypergeometric over measured genes; candidate genes p <= 0.01; GO Biological Process sizes 5-500',
      'transfer':'Spearman rank concordance across shared genes; AUC predicts direction among external genes with p <= 0.01',
      'caution':'Exploratory reanalysis; not a substitute for raw CEL-file normalization and preregistered confirmatory validation.'},
      'gene_level_summary':tissue_summary,'variance_priors':priors,
      'go_summary':{'tested_rows':len(go),'fdr_005':int((go.fdr<.05).sum()) if not go.empty else 0,'top':top_go},
      'transfer_summary':{'comparisons':len(transfer),'predictive':int((transfer.transfer_status=='predictive').sum()),
                          'rank_consistent':int((transfer.transfer_status=='rank-consistent').sum()),'top':best_transfer}}
    JSON_OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(payload,indent=2,ensure_ascii=False))


if __name__=='__main__': main()
