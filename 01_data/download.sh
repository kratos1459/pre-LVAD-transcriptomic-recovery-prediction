#!/usr/bin/env bash
# Stage 01: download public inputs (GEO series matrices + processed supplementary files).
# Usage: bash download.sh [DATA_DIR]   (default: ../data)
set -euo pipefail
D="${1:-$(dirname "$0")/../data}"; mkdir -p "$D"; cd "$D"
B=https://ftp.ncbi.nlm.nih.gov/geo/series
get() { # accession file-subdir filename
  local g=$1 sub=$2 f=$3 p=${1:0:${#1}-3}nnn
  [ -s "$f" ] || curl -sfL -o "$f" "$B/$p/$g/$sub/$f"
}
for g in GSE57338 GSE21610 GSE46224 GSE116250 GSE141910 GSE5406 GSE328059 GSE226314; do
  get $g matrix ${g}_series_matrix.txt.gz
done
get GSE46224  suppl GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz
get GSE116250 suppl GSE116250_rpkm.txt.gz
get GSE328059 suppl GSE328059_Derivation_Cohort_Repository_Analysis_R-Log_Values.csv.gz
get GSE328059 suppl GSE328059_Validation_Cohort_Repository_Analysis_R-Log_Values.csv.gz
# GSE57338 / GSE21610 / GSE141910 / GSE5406: expression is in the series matrix (or RAW.tar, see 01_data/README).
# Platform annotation (probe -> symbol) fetched by prepare_cohorts.py via GEOparse.
echo "done"; ls -lh
