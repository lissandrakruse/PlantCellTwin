#!/usr/bin/env python3
"""Targeted, study-aware validation for PlantCellTwin.

Uses NASA GeneLab processed differential-expression tables. GLDS-208 is a
root-zone reference dataset, not a spaceflight contrast. OSD-251 is an
in-flight fractional-gravity series. Results are intentionally panel-level
and retain the original study contrasts.
"""

from pathlib import Path
import json
import math
import pandas as pd
from scipy.stats import spearmanr, binomtest

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "nasa_expansion_data"
OUT = ROOT / "orbital-cell-twin"

panels = pd.read_csv(OUT / "expanded-marker-panels.csv")


def bh(values):
    s = pd.Series(values, dtype=float)
    valid = s.dropna().sort_values()
    n = len(valid)
    adjusted = pd.Series(float("nan"), index=s.index)
    if not n:
        return adjusted
    raw = valid.values * n / pd.Series(range(1, n + 1)).values
    raw = pd.Series(raw[::-1]).cummin()[::-1].clip(upper=1).values
    adjusted.loc[valid.index] = raw
    return adjusted


# GLDS-208: root zone I versus zone II, cross-platform confirmation.
r208 = pd.read_csv(DATA / "GLDS-208_rna_seq_differential_expression_GLbulkRNAseq.csv", low_memory=False)
a208 = pd.read_csv(DATA / "GLDS-208_array_differential_expression.csv", low_memory=False)
fc208 = "Log2fc_(ROOT ZONE I (0.5mm))v(ROOT ZONE II (1.5mm))"
q208 = "Adj.p.value_(ROOT ZONE I (0.5mm))v(ROOT ZONE II (1.5mm))"

r = panels.merge(r208[["TAIR", "SYMBOL", fc208, q208]], left_on="locus", right_on="TAIR", how="left")
r = r.rename(columns={fc208: "rna_log2fc_zoneI_vs_zoneII", q208: "rna_fdr"})
a = a208.groupby("TAIR", as_index=False).agg(
    array_log2fc_zoneI_vs_zoneII=(fc208, "median"), array_fdr=(q208, "min")
)
g208 = r.merge(a, left_on="locus", right_on="TAIR", how="left", suffixes=("", "_array"))
g208["rna_significant"] = g208["rna_fdr"] < 0.05
g208["array_significant"] = g208["array_fdr"] < 0.05
g208["platform_same_direction"] = (
    g208["rna_log2fc_zoneI_vs_zoneII"] * g208["array_log2fc_zoneI_vs_zoneII"] > 0
)
g208.to_csv(OUT / "glds208-root-zone-validation.csv", index=False)

# OSD-251: official uG versus 1G contrast plus an exploratory monotonic trend
# over group means at 0, 0.09, 0.18, 0.36, 0.57 and 1 G.
d251 = pd.read_csv(DATA / "GLDS-251_rna_seq_differential_expression_GLbulkRNAseq.csv", low_memory=False)
fc251 = "Log2fc_(Space Flight & uG)v(Space Flight & 1G by centrifugation)"
q251 = "Adj.p.value_(Space Flight & uG)v(Space Flight & 1G by centrifugation)"
gravities = [0.0, 0.09, 0.18, 0.36, 0.57, 1.0]
mean_cols = [
    "Group.Mean_(Space Flight & uG)",
    "Group.Mean_(Space Flight & 0.09G by centrifugation)",
    "Group.Mean_(Space Flight & 0.18G by centrifugation)",
    "Group.Mean_(Space Flight & 0.36G by centrifugation)",
    "Group.Mean_(Space Flight & 0.57G by centrifugation)",
    "Group.Mean_(Space Flight & 1G by centrifugation)",
]
x251 = panels.merge(d251[["TAIR", "SYMBOL", fc251, q251] + mean_cols], left_on="locus", right_on="TAIR", how="left")
x251 = x251.rename(columns={fc251: "log2fc_uG_vs_1G", q251: "fdr_uG_vs_1G"})
rhos, ps = [], []
for _, row in x251.iterrows():
    vals = [row[c] for c in mean_cols]
    if any(pd.isna(vals)):
        rhos.append(float("nan")); ps.append(float("nan"))
    else:
        rho, p = spearmanr(gravities, vals)
        rhos.append(rho); ps.append(p)
x251["gravity_spearman_rho"] = rhos
x251["gravity_trend_p"] = ps
x251["gravity_trend_fdr_panel"] = x251.groupby("panel")["gravity_trend_p"].transform(bh)
x251["same_direction_as_glds7_flight"] = x251["log2fc_uG_vs_1G"] * x251["glds7_root_log2fc"] > 0
x251.to_csv(OUT / "osd251-gravity-dose-validation.csv", index=False)

def finite(v):
    return None if pd.isna(v) or not math.isfinite(float(v)) else float(v)

summaries = []
for panel_name, sub in g208.groupby("panel"):
    paired = sub.dropna(subset=["rna_log2fc_zoneI_vs_zoneII", "array_log2fc_zoneI_vs_zoneII"])
    rho, p = spearmanr(paired["rna_log2fc_zoneI_vs_zoneII"], paired["array_log2fc_zoneI_vs_zoneII"])
    summaries.append({
        "study": "GLDS-208", "panel": panel_name, "role": "root-zone biological context",
        "n_expected": int(len(sub)), "n_available_rna": int(sub["rna_log2fc_zoneI_vs_zoneII"].notna().sum()),
        "n_rna_fdr_lt_0_05": int(sub["rna_significant"].sum()),
        "n_array_fdr_lt_0_05": int(sub["array_significant"].sum()),
        "n_cross_platform_same_direction": int(paired["platform_same_direction"].sum()),
        "cross_platform_spearman_rho": finite(rho), "cross_platform_spearman_p": finite(p),
    })

for panel_name, sub in x251.groupby("panel"):
    avail = sub.dropna(subset=["log2fc_uG_vs_1G"])
    pos_trend = int((sub["gravity_spearman_rho"] > 0).sum())
    sign_p = binomtest(pos_trend, int(sub["gravity_spearman_rho"].notna().sum()), 0.5).pvalue
    same = int(avail["same_direction_as_glds7_flight"].sum())
    summaries.append({
        "study": "OSD-251", "panel": panel_name, "role": "in-flight fractional-gravity response",
        "n_expected": int(len(sub)), "n_available": int(len(avail)),
        "n_uG_vs_1G_fdr_lt_0_05": int((avail["fdr_uG_vs_1G"] < 0.05).sum()),
        "n_uG_vs_1G_positive": int((avail["log2fc_uG_vs_1G"] > 0).sum()),
        "n_same_direction_as_glds7_flight": same,
        "n_positive_gravity_trends": pos_trend,
        "positive_trend_sign_test_p": float(sign_p),
        "median_gravity_spearman_rho": finite(sub["gravity_spearman_rho"].median()),
        "n_gene_trends_nominal_p_lt_0_05": int((sub["gravity_trend_p"] < 0.05).sum()),
        "n_gene_trends_panel_fdr_lt_0_05": int((sub["gravity_trend_fdr_panel"] < 0.05).sum()),
    })

summary_df = pd.DataFrame(summaries)
summary_df.to_csv(OUT / "nasa-expanded-validation-summary.csv", index=False)

payload = {
    "status": "completed numerical validation",
    "generated_from": [
        "NASA GeneLab GLDS-208 processed RNA-seq and microarray differential-expression tables",
        "NASA GeneLab OSD-251 processed RNA-seq differential-expression table",
    ],
    "claim_boundary": {
        "GLDS-208": "Root-zone reference only; it is not a spaceflight exposure comparison.",
        "OSD-251": "In-flight gravity-dose evidence in Ler-0 seedlings; it differs from GLDS-7 in ecotype, tissue composition, illumination, hardware and contrast.",
        "statistics": "Per-gene NASA FDR values are primary. Across-gene sign tests are exploratory because genes are not independent.",
    },
    "summary": summaries,
}
(OUT / "nasa-expanded-validation.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(summary_df.to_string(index=False))
