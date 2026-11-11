#!/usr/bin/env python3
"""Stage 04: novelty (PubMed via NCBI E-utilities), Enrichr over-representation, DGIdb approved-drug interactions.
Needs internet. Top 150 consensus genes (FDR<1e-6, concordant in >=7/8 cohorts; here >= n_cohorts-1)."""
import sys, time, json
sys.path.insert(0, "../02_meta_lasso_shap")
import requests, pandas as pd
from common import RES

res = pd.read_csv(RES / "meta_consensus.csv")
ncoh = res.k.max()
top = res[(res.fdr < 1e-6) & (res.n_concordant >= ncoh - 1)].head(150)
top.to_csv(RES / "top150_consensus.csv", index=False)

# --- PubMed counts (<=10 papers = under-studied)
Q = '({g}[Title/Abstract]) AND ("heart failure"[Title/Abstract] OR cardiomyopathy[Title/Abstract] OR "ventricular remodeling"[Title/Abstract] OR LVAD[Title/Abstract] OR "left ventricular assist device"[Title/Abstract])'
cnt = {}
for g in top.gene:
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={"db": "pubmed", "term": Q.format(g=g), "retmode": "json"}, timeout=30)
    cnt[g] = int(r.json()["esearchresult"]["count"])
    time.sleep(0.4)
top = top.assign(pubmed_hf_papers=top.gene.map(cnt), understudied=lambda d: d.pubmed_hf_papers <= 10)
print("under-studied:", int(top.understudied.sum()), "of", len(top))

# --- Enrichr
LIBS = ["GO_Biological_Process_2023", "Reactome_2022", "KEGG_2021_Human", "MSigDB_Hallmark_2020"]
import gseapy as gp
rows = []
for name, genes in [("up", top[top.g > 0].gene.tolist()), ("down", top[top.g < 0].gene.tolist())]:
    e = gp.enrichr(gene_list=genes, gene_sets=LIBS, organism="human", outdir=None).results
    e["direction"] = name
    rows.append(e)
pd.concat(rows).to_csv(RES / "enrichr.csv", index=False)

# --- DGIdb (GraphQL v5)
q = """query($g:[String!]){genes(names:$g){nodes{name interactions{drug{name approved} interactionTypes{type directionality} }}}}"""
out = []
for i in range(0, len(top), 50):
    r = requests.post("https://dgidb.org/api/graphql", json={"query": q, "variables": {"g": top.gene.tolist()[i:i + 50]}}, timeout=60).json()
    for n in r["data"]["genes"]["nodes"]:
        for it in n["interactions"]:
            if it["drug"]["approved"]:
                out.append(dict(gene=n["name"], drug=it["drug"]["name"], types=";".join(t["type"] or "" for t in it["interactionTypes"])))
pd.DataFrame(out).to_csv(RES / "dgidb_approved.csv", index=False)
top.to_csv(RES / "top150_annotated.csv", index=False)
