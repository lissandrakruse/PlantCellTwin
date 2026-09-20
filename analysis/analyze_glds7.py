import gzip
import json
import numpy as np
import pandas as pd
from scipy import stats
from config import DATA_DIR, RESULTS_DIR, require


MATRIX = DATA_DIR / 'GSE56659_series_matrix.txt.gz'
ANNOT = DATA_DIR / 'GPL198.annot.gz'
OUT = RESULTS_DIR / 'glds7-results.json'
CSV = RESULTS_DIR / 'glds7-top-genes.csv'


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def read_matrix():
    require(MATRIX, 'GEO series matrix')
    titles = None
    accessions = None
    header = None
    rows = []
    in_table = False
    with gzip.open(MATRIX, 'rt', errors='replace') as fh:
        for line in fh:
            if line.startswith('!Sample_title'):
                titles = [x.strip('"') for x in line.rstrip().split('\t')[1:]]
            elif line.startswith('!Sample_geo_accession'):
                accessions = [x.strip('"') for x in line.rstrip().split('\t')[1:]]
            elif line.startswith('!series_matrix_table_begin'):
                in_table = True
            elif line.startswith('!series_matrix_table_end'):
                break
            elif in_table and header is None:
                header = [x.strip('"') for x in line.rstrip().split('\t')]
            elif in_table:
                parts = line.rstrip().split('\t')
                rows.append([parts[0].strip('"')] + [float(x) for x in parts[1:]])
    frame = pd.DataFrame(rows, columns=header).set_index('ID_REF')
    meta = pd.DataFrame({'accession': accessions, 'title': titles})
    meta['condition'] = np.where(meta.title.str.startswith('FLIGHT'), 'flight', 'ground')
    meta['tissue'] = 'whole plant'
    meta.loc[meta.title.str.contains('SHOOTS'), 'tissue'] = 'shoot'
    meta.loc[meta.title.str.contains('HYPOCOTYLS'), 'tissue'] = 'hypocotyl'
    meta.loc[meta.title.str.contains('ROOTS'), 'tissue'] = 'root'
    return frame, meta


def read_annotation():
    require(ANNOT, 'GEO platform annotation')
    with gzip.open(ANNOT, 'rt', errors='replace') as fh:
        for line in fh:
            if line.startswith('!platform_table_begin'):
                break
        return pd.read_csv(fh, sep='\t', dtype=str)


def clean(value):
    if pd.isna(value):
        return ''
    return str(value).split(' /// ')[0].strip()


def analyze():
    expr, meta = read_matrix()
    ann = read_annotation().set_index('ID')
    symbols = ann.get('Gene symbol', pd.Series(index=ann.index, dtype=str))
    titles = ann.get('Gene title', pd.Series(index=ann.index, dtype=str))
    loci = ann.get('Platform_ORF', pd.Series(index=ann.index, dtype=str))
    log_expr = np.log2(expr + 1)
    results = {'dataset': {
        'n_samples': int(meta.shape[0]), 'n_probes': int(expr.shape[0]),
        'platform': 'Affymetrix Arabidopsis ATH1-121501 (GPL198)',
        'contrast': 'spaceflight vs matched 1 g ground control',
        'age': '12 days after light-stimulated germination',
        'mission': 'STS-131 / ISS, TAGES experiment',
        'accessions': ['NASA GLDS-7', 'GEO GSE56659']
    }, 'tissues': {}}
    all_rows = []
    for tissue in ['whole plant', 'shoot', 'hypocotyl', 'root']:
        ids_f = meta.loc[(meta.tissue == tissue) & (meta.condition == 'flight'), 'accession'].tolist()
        ids_g = meta.loc[(meta.tissue == tissue) & (meta.condition == 'ground'), 'accession'].tolist()
        f = log_expr[ids_f].to_numpy()
        g = log_expr[ids_g].to_numpy()
        lfc = np.nanmean(f, axis=1) - np.nanmean(g, axis=1)
        _, p = stats.ttest_ind(f, g, axis=1, equal_var=False, nan_policy='omit')
        p = np.nan_to_num(p, nan=1.0)
        fdr = bh(p)
        table = pd.DataFrame({'probe': expr.index, 'log2fc': lfc, 'p': p, 'fdr': fdr})
        table['symbol'] = [clean(symbols.get(i, '')) for i in table.probe]
        table['title'] = [clean(titles.get(i, '')) for i in table.probe]
        table['locus'] = [clean(loci.get(i, '')).upper() for i in table.probe]
        table['abs_fc'] = table.log2fc.abs()
        sig = table[(table.fdr < 0.05) & (table.abs_fc >= 1)]
        candidates = table[(table.p < 0.05) & (table.abs_fc >= 1)]
        ranked = table.sort_values(['fdr', 'abs_fc'], ascending=[True, False]).head(12)
        top = []
        for _, row in ranked.iterrows():
            top.append({'probe': row.probe, 'locus': row.locus, 'symbol': row.symbol or row.probe,
                        'title': row.title, 'log2fc': round(float(row.log2fc), 3),
                        'fdr': round(float(row.fdr), 6)})
            all_rows.append({'tissue': tissue, **top[-1]})
        results['tissues'][tissue] = {
            'n_flight': len(ids_f), 'n_ground': len(ids_g),
            'significant_total': int(sig.shape[0]),
            'up': int((sig.log2fc > 0).sum()), 'down': int((sig.log2fc < 0).sum()),
            'exploratory_candidates': int(candidates.shape[0]),
            'median_abs_log2fc': round(float(np.median(np.abs(lfc))), 3),
            'top_genes': top
        }
    OUT.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding='utf-8')
    pd.DataFrame(all_rows).to_csv(CSV, index=False)


if __name__ == '__main__':
    analyze()
