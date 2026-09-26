#!/usr/bin/env bash
# End-to-end validation run, first spike the deposited truth sets into clean synthetic
# BAMs then calls MNVista, lastly score the calls against the truth.
#
# Usage:
#   run_validation.sh <work_dir> <uppercased_ref.fa>
#
# Expects <work_dir>/bams/S1.bam ... S5.bam from make_synthetic_bams.py, and
# $BAMEDIT / $MNVISTA pointing at the two binaries
#
# Resumable: work already done is skipped, and each spiked BAM is deleted once
# it has been called, so the run never needs space for all 15
#
# BY_REGION=1 calls MNVista once per target region instead of once per sample.
# MNVista processes each window independently, so the two give the same result.
# Per-region calls finishes in ~ under a minute, for machines that
# cannot hold a long job or are less powerful.

set -euo pipefail

WORK=${1:?work dir}
REF=${2:?uppercased reference fasta}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE="$(dirname "$HERE")"
BAMEDIT=${BAMEDIT:-bamedit}
MNVISTA=${MNVISTA:-mnvista}
BY_REGION=${BY_REGION:-1}
JOBS=${JOBS:-1}
THREADS=${THREADS:-4}
REGIONS="$BUNDLE/spike_in_design/regions_se.bed"

# MNVista settings
MNV_ARGS=(--read-length 150 --read-quality 30 --max-mnv-size -1
          --bayes-prior-mnv 0.5 --bayes-p-error 0.5 --min-bayesian 0.95
          --min-phi 0 --min-jaccard 0 --min-vaf-mnv 0.000001
          --min-vrd-snv 1 --min-vrd-mnv 1 --min-vaf-snv 0 --max-vaf-snv 1)

mkdir -p "$WORK/spiked" "$WORK/mnvista_out" "$WORK/parts"

split_vcf() {   #If per-region is enables then split the VCF
  python3 - "$WORK" "$1" "$REGIONS" <<'PY'
import os, sys
work, name, regions = sys.argv[1], sys.argv[2], sys.argv[3]
vcf = f"{work}/spiked/{name}.vcf"
out = f"{work}/parts/{name}"
os.makedirs(out, exist_ok=True)
lines = open(vcf).readlines()
hdr = [l for l in lines if l.startswith("#")]
rows = [l for l in lines if not l.startswith("#")]
for i, r in enumerate(l.split() for l in open(regions)):
    c, s, e = r[0], int(r[1]), int(r[2])
    sub = [l for l in rows
           if l.split("\t")[0] == c and s <= int(l.split("\t")[1]) < e]
    open(f"{out}/r{i:02d}.vcf", "w").writelines(hdr + sub)
PY
}

for S in S1 S2 S3 S4 S5; do
  for A in noise0 noise1 noise10; do
    N="${S}_${A}"
    [[ -f "$WORK/mnvista_out/$N.ok" ]] && continue

    if [[ -f "$WORK/spiked/$N.bam" ]] && ! python3 -c "
import pysam,sys
pysam.AlignmentFile(sys.argv[1]).close()" "$WORK/spiked/$N.bam" 2>/dev/null; then
      echo "[spike] $N: previous output truncated, redoing"
      rm -f "$WORK/spiked/$N.bam" "$WORK/spiked/$N.bam.bai" "$WORK/spiked/$N.vcf"
    fi

    if [[ ! -f "$WORK/spiked/$N.bam" || ! -f "$WORK/spiked/$N.vcf" ]]; then
      echo "[spike] $N"
      rm -f "$WORK/spiked/$N.bam" "$WORK/spiked/$N.bam.bai" "$WORK/spiked/$N.vcf"
      # -S 1 is BamEdit's default seed read selection is deterministic.
      "$BAMEDIT" "$BUNDLE/truth_sets/beds/$N.bed" "$REF" "$WORK/bams/$S.bam" \
                 "$WORK/spiked" -O "$N" -S 1 -T 4
    fi
    [[ -f "$WORK/spiked/$N.bam.bai" ]] || \
      python3 -c "import pysam,sys;pysam.index(sys.argv[1])" "$WORK/spiked/$N.bam"

    if [[ "$BY_REGION" == "1" ]]; then
      [[ -d "$WORK/parts/$N" ]] || split_vcf "$N"
      # Regions are independent, so run $JOBS of them 
      for V in "$WORK/parts/$N"/r*.vcf; do
        R=$(basename "$V" .vcf)
        [[ -f "$WORK/parts/$N/$R.csv" ]] && continue
        echo "[call] $N $R"
        "$MNVISTA" "$WORK/spiked/$N.bam" "$V" "$WORK/parts/$N" -O "$R" \
                   "${MNV_ARGS[@]}" --threads "$THREADS" &
        while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n; done
      done
      wait
      for V in "$WORK/parts/$N"/r*.vcf; do
        R=$(basename "$V" .vcf)
        grep -q "Finished!" "$WORK/parts/$N/$R.log"
      done
      python3 - "$WORK/parts/$N" "$WORK/mnvista_out/$N.csv" <<'PY'
import glob, os, sys
src, dst = sys.argv[1], sys.argv[2]
parts = sorted(p for p in glob.glob(os.path.join(src, "r*.csv"))
               if "filtered" not in p)
with open(dst, "w") as out:
    for i, p in enumerate(parts):
        with open(p) as fh:
            head = fh.readline()
            if i == 0:
                out.write(head)
            out.writelines(fh)
PY
    else
      echo "[call] $N"
      "$MNVISTA" "$WORK/spiked/$N.bam" "$WORK/spiked/$N.vcf" \
                 "$WORK/mnvista_out" -O "$N" "${MNV_ARGS[@]}" \
                 --threads "$THREADS"
      grep -q "Finished!" "$WORK/mnvista_out/$N.log"
    fi

    touch "$WORK/mnvista_out/$N.ok"
    rm -f "$WORK/spiked/$N.bam" "$WORK/spiked/$N.bam.bai"
  done
done

python3 "$HERE/score_calls.py" "$BUNDLE/truth_sets/mnv" \
        "$WORK/mnvista_out" "$WORK/RESULTS.csv"
