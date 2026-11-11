"""Ensembl gene id -> HGNC-style symbol, built from Ensembl+symbol tables already shipped with
GSE116250 (Common_name) and GSE328059 (Associated Gene Name / gene_name). No external download."""
from pathlib import Path
import pandas as pd
D = Path(__file__).resolve().parent.parent / "data"
_m = None


def _build():
    a = pd.read_csv(D / "GSE116250_rpkm.txt.gz", sep="\t", usecols=[0, 1])
    a.columns = ["id", "sym"]
    b = pd.read_csv(D / "GSE328059_Derivation_Cohort_Repository_Analysis_R-Log_Values.csv.gz", usecols=[0, 1])
    b.columns = ["id", "sym"]
    c = pd.read_csv(D / "GSE328059_Validation_Cohort_Repository_Analysis_R-Log_Values.csv.gz", usecols=[0, 1])
    c.columns = ["id", "sym"]
    t = pd.concat([a, b, c]).dropna().drop_duplicates("id")
    t["id"] = t["id"].str.split(".").str[0]
    return t.set_index("id")["sym"]


def ensembl_to_symbol(ids):
    global _m
    if _m is None:
        _m = _build()
    ids = pd.Index(ids).astype(str).str.split(".").str[0]
    return _m.reindex(ids).values
