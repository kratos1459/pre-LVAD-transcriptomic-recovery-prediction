#!/usr/bin/env bash
# Stage 06a: PRJNA1312133 (ENA) -> Salmon against GENCODE v46. Streaming: FASTQ deleted after each run.
# Only pre-LVAD HF (*_HF_*) and non-failing (*_NF) runs. ~9 GB/run; P parallel runs
# (the link is capped ~5 MB/s in total, so expect ~10-12 h).
# Usage: bash quantify_prjna1312133.sh WORKDIR [THREADS_PER_JOB=2] [PARALLEL=8]   (salmon on PATH, curl)
set -euo pipefail
W=${1:?workdir}; T=${2:-2}; P=${3:-8}; mkdir -p "$W"/{ref,fq,quant,logs}; cd "$W"
[ -s ref/gencode.v46.transcripts.fa.gz ] || curl -sfL -o ref/gencode.v46.transcripts.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz
[ -d ref/idx ] || salmon index -t ref/gencode.v46.transcripts.fa.gz -i ref/idx -p 12 --gencode
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1312133&result=read_run&fields=run_accession,sample_title,fastq_ftp&format=tsv" | tail -n +2 > runs.tsv
grep -E "_HF_|_NF" runs.tsv > todo.tsv
cat > one.sh <<'EOF'
#!/usr/bin/env bash
run=$1; title=$2; ftp=$3; T=$4
[ -s quant/$title/quant.sf ] && exit 0
IFS=';' read -r f1 f2 <<< "$ftp"
curl -sfL --retry 8 -C - -o fq/${run}_1.fastq.gz "https://$f1"
curl -sfL --retry 8 -C - -o fq/${run}_2.fastq.gz "https://$f2"
if gzip -t fq/${run}_1.fastq.gz && gzip -t fq/${run}_2.fastq.gz && salmon quant -i ref/idx -l A -1 fq/${run}_1.fastq.gz -2 fq/${run}_2.fastq.gz -p "$T" --validateMappings --gcBias -o quant/$title > logs/$title.log 2>&1 && [ -s quant/$title/quant.sf ]; then echo "done $title"; else echo "FAILED $title"; fi
rm -f fq/${run}_*.fastq.gz
EOF
chmod +x one.sh
cut -f1-3 todo.tsv | tr '\t' '\n' | paste - - - | xargs -P "$P" -L1 bash -c './one.sh "$0" "$1" "$2" '"$T"
echo ALL_DONE
