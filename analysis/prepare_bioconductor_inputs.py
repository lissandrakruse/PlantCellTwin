import re
import pandas as pd

from analyze_glds7 import read_matrix, read_annotation, clean
from config import DATA_DIR

ROOT = DATA_DIR

_, meta = read_matrix()
meta['plate'] = meta['title'].str.extract(r'plate(\d)', expand=False).astype(int)
meta.to_csv(ROOT / 'glds7_sample_metadata.csv', index=False)

ann = read_annotation()
mapping = pd.DataFrame({
    'probe': ann['ID'],
    'locus': ann['Platform_ORF'].map(clean).str.upper(),
    'symbol': ann.get('Gene symbol', '').map(clean),
    'title': ann.get('Gene title', '').map(clean),
})
mapping = mapping[mapping['locus'].str.match(r'^AT[1-5CM]G\d{5}$', na=False)]
mapping.to_csv(ROOT / 'glds7_probe_locus_map.csv', index=False)
