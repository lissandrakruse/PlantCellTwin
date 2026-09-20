import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from analyze_glds7 import read_matrix, read_annotation, clean, bh

ROOT = Path('/workspace/scratch/1b7c2e5453f2')
BOOK = ROOT / 'arabidopsis_MOESM5.xlsx'
OUT = ROOT / 'orbital-cell-twin' / 'nasa-cross-validation.json'
CSV = ROOT / 'orbital-cell-twin' / 'nasa-cross-validation.csv'


SHEETS = {
    'GLDS-17': ('Arabidopsis Gene stable ID', 'Log2fc_(Space Flight & Seedlings)v(Ground Control & Seedlings)', 'P.value_(Space Flight & Seedlings)v(Ground Control & Seedlings)'),
    'GLDS-37': ('Arabidopsis Gene stable ID', 'GLDS-37_Log2fc_(Col-0_FLT)v(Col-0_GC)', 'GLDS-37_P.value_(Col-0_FLT)v(Col-0_GC)'),
    'GLDS-38': ('Arabidopsis Gene stable ID.1', 'GLDS-38_Log2fc_(FLT)v(GC)', 'GLDS-38_P.value_(FLT)v(GC)'),
    'GLDS-44': ('Arabidopsis Gene stable ID', 'GLDS44_Log2fc_(wild type & Space Flight)v(wild type & Ground Control)', 'GLDS44_P.value_(wild type & Space Flight)v(wild type & Ground Control)'),
    'GLDS-121': ('Arabidopsis Gene stable ID', 'GLDS-121 Log2fc_(Space Flight)v(Ground Control)', 'GLDS-121 P.value_(Space Flight)v(Ground Control)'),
}


def external_sets():
    out = {}
    for sheet, (id_col, fc_col, p_col) in SHEETS.items():
        df = pd.read_excel(BOOK, sheet_name=sheet, header=1)
        ids = df[id_col]
        if isinstance(ids, pd.DataFrame):
            ids = ids.iloc[:, -1]
        part = pd.DataFrame({'locus': ids, 'log2fc': pd.to_numeric(df[fc_col], errors='coerce'), 'p': pd.to_numeric(df[p_col], errors='coerce')})
        part['locus'] = part.locus.astype(str).str.upper().str.strip()
        part = part[part.locus.str.match(r'^AT[1-5CM]G\d{5}$', na=False)].dropna(subset=['log2fc', 'p'])
        part = part.sort_values('p').drop_duplicates('locus')
        out[sheet] = part.set_index('locus')
    return out


def glds7_stats():
    expr, meta = read_matrix()
    ann = read_annotation().set_index('ID')
    loci = ann.get('Platform_ORF', pd.Series(index=ann.index, dtype=str)).map(clean).str.upper()
    log_expr = np.log2(expr + 1)
    out = {}
    for tissue in ['whole plant', 'shoot', 'hypocotyl', 'root']:
        ids_f = meta.loc[(meta.tissue == tissue) & (meta.condition == 'flight'), 'accession'].tolist()
        ids_g = meta.loc[(meta.tissue == tissue) & (meta.condition == 'ground'), 'accession'].tolist()
        f, g = log_expr[ids_f].to_numpy(), log_expr[ids_g].to_numpy()
        lfc = np.nanmean(f, axis=1) - np.nanmean(g, axis=1)
        _, p = stats.ttest_ind(f, g, axis=1, equal_var=False, nan_policy='omit')
        frame = pd.DataFrame({'locus': loci.reindex(expr.index).values, 'log2fc': lfc, 'p': np.nan_to_num(p, nan=1.0)})
        frame = frame[frame.locus.str.match(r'^AT[1-5CM]G\d{5}$', na=False)].sort_values('p').drop_duplicates('locus')
        out[tissue] = frame.set_index('locus')
    return out


def main():
    g7, ext = glds7_stats(), external_sets()
    universe = len(set().union(*[set(x.index) for x in g7.values()], *[set(x.index) for x in ext.values()]))
    rows = []
    for tissue, a_all in g7.items():
        a = a_all[a_all.p <= 0.01]
        for study, b_all in ext.items():
            b = b_all[b_all.p <= 0.01]
            common = sorted(set(a.index) & set(b.index))
            same = sum(np.sign(a.loc[g, 'log2fc']) == np.sign(b.loc[g, 'log2fc']) for g in common)
            direction_rate = same / len(common) if common else np.nan
            overlap_p = stats.hypergeom.sf(len(common)-1, universe, len(a), len(b)) if common else 1.0
            direction_p = stats.binomtest(same, len(common), .5, alternative='greater').pvalue if common else 1.0
            rows.append({'tissue': tissue, 'study': study, 'glds7_candidates': len(a), 'external_candidates': len(b),
                         'overlap': len(common), 'direction_agreement': direction_rate,
                         'overlap_p': overlap_p, 'direction_p': direction_p, 'genes': common[:20]})
    overlap_fdr = bh([r['overlap_p'] for r in rows])
    direction_fdr = bh([r['direction_p'] for r in rows])
    for r, op, dp in zip(rows, overlap_fdr, direction_fdr):
        r['overlap_fdr'] = float(op); r['direction_fdr'] = float(dp)
        r['status'] = 'replicated' if op < .05 and dp < .05 and r['direction_agreement'] >= .6 else ('partial' if r['overlap'] >= 5 and r['direction_agreement'] >= .5 else 'not replicated')
    best = sorted(rows, key=lambda r: (r['status'] != 'replicated', r['status'] != 'partial', r['overlap_fdr'], -r['overlap']))
    payload = {
        'method': {'threshold': 'p ≤ 0.01 in each study', 'multiple_testing': 'BH across 20 tissue–study comparisons',
                   'overlap_test': 'one-sided hypergeometric', 'direction_test': 'one-sided exact binomial',
                   'source': 'Barker et al. 2023 Supplementary Data 4; GeneLab common-pipeline columns',
                   'universe_genes': universe},
        'summary': {'comparisons': len(rows), 'replicated': sum(r['status']=='replicated' for r in rows),
                    'partial': sum(r['status']=='partial' for r in rows), 'not_replicated': sum(r['status']=='not replicated' for r in rows)},
        'comparisons': best
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')
    flat = [{k: (','.join(v) if isinstance(v, list) else v) for k,v in r.items()} for r in best]
    pd.DataFrame(flat).to_csv(CSV, index=False)


if __name__ == '__main__':
    main()
