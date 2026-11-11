#!/usr/bin/env python3
"""Stage 01: harmonise cohorts to gene-symbol x sample log2-scale matrices + sample metadata.
Output (../results/cohorts/): <cohort>_expr.csv.gz, all_meta.csv   (HF = failing/pre-LVAD, NF = donor)
Rules (manuscript Methods): probes collapsed to max-mean symbol; RNA-seq log2-transformed if max>100.
"""
import gzip, io, re, tarfile, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
D, OUT = ROOT / "data", ROOT / "results" / "cohorts"
OUT.mkdir(parents=True, exist_ok=True)


def series_meta(g):
    rows = {}
    for l in gzip.open(D / f"{g}_series_matrix.txt.gz", "rt"):
        if l.startswith("!Sample_") and not l.startswith("!Sample_supplementary"):
            k = l.split("\t")[0][1:]
            v = [x.strip().strip('"') for x in l.rstrip("\n").split("\t")[1:]]
            i = 1
            while f"{k}#{i}" in rows:
                i += 1
            rows[k if k not in rows else f"{k}#{i}"] = v
    return pd.DataFrame(rows)


def series_expr(g):
    lines = []
    on = False
    for l in gzip.open(D / f"{g}_series_matrix.txt.gz", "rt"):
        if l.startswith("!series_matrix_table_begin"):
            on = True
            continue
        if l.startswith("!series_matrix_table_end"):
            break
        if on:
            lines.append(l)
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0).rename(columns=lambda c: c.strip('"')).rename(index=lambda i: str(i).strip('"'))


def probe_map(gpl):
    f = D / "plat" / f"{gpl}.annot.gz"
    txt = [l for l in gzip.open(f, "rt", errors="ignore") if not l.startswith(("!", "^", "#"))]
    df = pd.read_csv(io.StringIO("".join(txt)), sep="\t", low_memory=False)
    col = [c for c in df.columns if c.lower() in ("gene symbol", "gene_symbol")][0]
    s = df[col].astype(str).str.split("///").str[0].str.strip()
    s = s.where(~s.isin(["nan", "", "---"]))
    return pd.Series(s.values, index=df["ID"].astype(str).values).dropna()


def collapse(expr, mp):
    """probe -> symbol, keep max-mean probe per symbol."""
    expr = expr.copy()
    expr.index = expr.index.astype(str)
    sym = mp.reindex(expr.index)
    ok = sym.notna()
    e = expr[ok]
    e["sym"] = sym[ok].values
    e["m"] = e.drop(columns="sym").mean(axis=1)
    e = e.sort_values("m", ascending=False).drop_duplicates("sym").drop(columns="m").set_index("sym")
    return e


def log2_if_needed(m):
    return np.log2(m.clip(lower=0) + 1) if np.nanmax(m.values) > 100 else m


def collapse_symbol_rows(df):
    df = df[df.index.notna() & (df.index != "")]
    mm = df.mean(axis=1)
    df = df.assign(_m=mm).sort_values("_m", ascending=False)
    return df[~df.index.duplicated()].drop(columns="_m")


cohorts, metas = {}, []


def add(name, expr, meta):
    meta = meta.copy()
    meta["cohort"] = name
    meta.index = meta["sample"]
    expr = collapse_symbol_rows(expr)
    expr = expr[meta["sample"]]
    expr.to_csv(OUT / f"{name}_expr.csv.gz")
    metas.append(meta)
    print(name, expr.shape, meta["group"].value_counts().to_dict())


# ---- GSE57338 (Affymetrix HuGene 1.1 ST, GPL11532) ----
m = series_meta("GSE57338")
hf = m["Sample_characteristics_ch1#1"].str.contains("yes") if "Sample_characteristics_ch1#1" in m else None
hfcol = [c for c in m.columns if m[c].astype(str).str.contains("heart failure:").any()][0]
dis = [c for c in m.columns if m[c].astype(str).str.contains("disease status:").any()][0]
e = collapse(series_expr("GSE57338"), probe_map("GPL11532"))
meta = pd.DataFrame({"sample": m["Sample_geo_accession"], "group": np.where(m[hfcol].str.contains("yes"), "HF", "NF"),
                     "etiology": m[dis].str.replace("disease status: ", "")})
add("GSE57338", e, meta)

# ---- GSE21610 (GPL570): NF vs pre-VAD HF; post-VAD kept separately for reversal ----
m = series_meta("GSE21610")
vcol = [c for c in m.columns if m[c].astype(str).str.contains("vad status:").any()][0]
dcol = [c for c in m.columns if m[c].astype(str).str.contains("disease status:").any()][0]
v = m[vcol].str.replace("vad status: ", "")
e = collapse(series_expr("GSE21610"), probe_map("GPL570"))
grp = v.map({"unsupported": "NF", "before VAD support": "HF", "after VAD support": "POST"})
meta = pd.DataFrame({"sample": m["Sample_geo_accession"], "group": grp, "etiology": m[dcol].str.replace("disease status: ", "")})
add("GSE21610", e, meta)

# ---- GSE5406 (GPL96 = subset of GPL570 probe IDs) ----
m = series_meta("GSE5406")
ch = m["Sample_characteristics_ch1"]
e = collapse(series_expr("GSE5406"), probe_map("GPL570"))
meta = pd.DataFrame({"sample": m["Sample_geo_accession"], "group": np.where(ch.str.contains("donor"), "NF", "HF"),
                     "etiology": np.where(ch.str.contains("ischemic cardiomyopathy"), "ischemic", np.where(ch.str.contains("donor"), "none", "idiopathic DCM"))})
add("GSE5406", e, meta)

# ---- GSE46224 (RNA-seq RPKM; columns ICM, ICM+LVAD(post), NF, NICM, NICM+LVAD) ----
r = pd.read_csv(D / "GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz", sep="\t", index_col=0)
r = log2_if_needed(r)
cols = list(r.columns)
grp = ["POST" if "LVAD" in c else ("NF" if c.startswith("NF") else "HF") for c in cols]
et = ["ischemic" if c.startswith("ICM") else ("none" if c.startswith("NF") else "non-ischemic") for c in cols]
add("GSE46224", r, pd.DataFrame({"sample": cols, "group": grp, "etiology": et}))

# ---- GSE116250 (RNA-seq, log2 RPKM-like; Ensembl + common name) ----
r = pd.read_csv(D / "GSE116250_rpkm.txt.gz", sep="\t", index_col=0)
sym = r.pop("Common_name")
r.index = sym.values
r = log2_if_needed(r)
cols = list(r.columns)
grp = ["NF" if c.startswith("NF") else "HF" for c in cols]
et = ["none" if c.startswith("NF") else ("ischemic" if c.startswith("ICM") else "DCM") for c in cols]
add("GSE116250", r, pd.DataFrame({"sample": cols, "group": grp, "etiology": et}))

# ---- GSE141910 (per-sample CSVs, Ensembl ids; DCM/HCM/PPCM vs donor) ----
m = series_meta("GSE141910")
ecol = [c for c in m.columns if m[c].astype(str).str.contains("etiology:").any()][0]
frames = {}
for f in sorted((D / "g141910").glob("*.csv.gz")):
    s = f.name.split("_")[1].replace(".csv.gz", "")
    frames[s] = pd.read_csv(f, index_col=0).iloc[:, 0]
x = pd.DataFrame(frames)
from ens_map import ensembl_to_symbol  # noqa: E402  (built by build_ens_map.py from GSE116250 table)
x.index = ensembl_to_symbol(x.index)
x = log2_if_needed(x)
et = m[ecol].str.replace("etiology: ", "")
meta = pd.DataFrame({"sample": m["Sample_title"], "group": np.where(et.str.contains("Non-Failing"), "NF", "HF"), "etiology": et})
add("GSE141910", x, meta)

# ---- GSE328059 (derivation + validation; labels live in GEO sample titles) ----
m = series_meta("GSE328059")
lib = m["Sample_description"].str.replace("Library name: ", "")
tit = m["Sample_title"]
lab = tit.str.extract(r"(Derivation|Validation) Cohort, (PRE-NR|PRE-R|POST-NR|POST-R|DONOR|Patient)")
dev = pd.read_csv(D / "GSE328059_Derivation_Cohort_Repository_Analysis_R-Log_Values.csv.gz").drop(columns="Gene").set_index("Associated Gene Name")
val = pd.read_csv(D / "GSE328059_Validation_Cohort_Repository_Analysis_R-Log_Values.csv.gz").drop(columns="ensembl_id")
val = val.loc[:, ~val.columns.str.startswith("Unnamed")].set_index("gene_name")
dev, val = collapse_symbol_rows(dev), collapse_symbol_rows(val)
both = dev.join(val, how="inner")
rows = []
for i in range(len(m)):
    rows.append(dict(sample=lib.iloc[i], gsm=m["Sample_geo_accession"].iloc[i], cohort_part=lab.iloc[i, 0], label=lab.iloc[i, 1],
                     batch=re.match(r"(\d+)X", lib.iloc[i]).group(1)))
meta = pd.DataFrame(rows)
meta["group"] = meta["label"].map({"DONOR": "NF", "PRE-R": "HF", "PRE-NR": "HF", "POST-R": "POST", "POST-NR": "POST", "Patient": "UNLABELLED"})
meta["recovery"] = meta["label"].map({"PRE-R": "R", "PRE-NR": "NR", "POST-R": "R", "POST-NR": "NR"})
meta["timepoint"] = meta["label"].map(lambda s: "pre" if str(s).startswith("PRE") else ("post" if str(s).startswith("POST") else ("donor" if s == "DONOR" else "unlabelled")))
add("GSE328059", both, meta)

pd.concat(metas, ignore_index=True).to_csv(OUT / "all_meta.csv", index=False)
print("saved", OUT)
