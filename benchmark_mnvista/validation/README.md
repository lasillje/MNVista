# Self-contained rerun of the benchmark

1. build five clean synthetic BAMs over the 10 target regions, containing no
   variants at all
2. spike the deposited BEDs into them with BamEdit
3. call MNVista on the spiked BAMs, using BamEdit's own VCF as the candidate
   list and the exact settings from the manuscript run
4. score the calls against the deposited MNV truth sets.

Nothing in this folder depend on patient data, or on any file that
is not in this repository. 

Note that the benchmark can take quite a while to complete in its entirety, as all sequence data has to be regenerated.

## Result

**8,008 of 8,008 spiked MNVs recalled.**

| | truth MNVs | calls | TP | FP | FN | precision | recall | F1 |
|---|---|---|---|---|---|---|---|---|
| All 15 runs | 8,008 | 8,049 | 8,008 | 41 | **0** | 0.9949 | **1.0000** | 0.9974 |

Per sample and arm, see `results/RESULTS.csv`. Recall is 1.000 in all fifteen.
Precision ranges from 0.974 (S2, 10% noise) to 1.000 (S2/S3/S4).

### False positives

Every false positive is a pairing error among the lowest-abundance spiked SNVs,
and every one is barely supported:

| Threshold | False positives remaining | True positives lost |
|---|---|---|
| ≥ 2 fragments | 11 / 41 | 0 / 8,008 |
| **≥ 3 fragments** | **1 / 41** | **0 / 8,008** |
| ≥ 4 fragments | 0 / 41 | 1 / 8,008 |


## Running the benchmark re-run

```bash
# 1. an uppercased reference fasta (no softmasks)
awk '/^>/ {print; next} {print toupper($0)}' hg38.fa > hg38.upper.fa
samtools faidx hg38.upper.fa

# 2. five clean BAMs at 18,000x over the 10 regions (~4M reads each)
python make_synthetic_bams.py hg38.upper.fa ../spike_in_design/regions_se.bed work/bams

# 3. spike, call, score
BAMEDIT=/path/to/bamedit MNVISTA=/path/to/mnvista \
  ./run_validation.sh work hg38.upper.fa

# 4. If you did not build the binaries in another location, i.e. MNVista at
#    MNVista/bin/MNVista and BamEdit at benchmark_mnvista/bamedit/bin/bamedit,
#    and you run this from benchmark_mnvista/validation/, you can use:

BAMEDIT=../bamedit/bin/bamedit MNVISTA=../../bin/MNVista \
  ./run_validation.sh work hg38.upper.fa
```

`run_validation.sh` is resumable and deletes each spiked BAM once it has been
called, so the run never needs space for all fifteen BAMs at once. `BY_REGION=1`
calls MNVista once per target region rather than once per sample, for machines
that cannot hold a long job (e.g. in a cluster environment, or outdated PC specs); `JOBS`/`THREADS` control resource concurrency if ran on a cluster.

| Script | Result |
|---|---|
| `make_synthetic_bams.py` | Generates clean BAMs. |
| `score_calls.py` | Scores MNVista output against `../truth_sets/mnv/`. |
| `run_validation.sh` | Spike, call, and score. |

## The synthetic BAMs

Reads are copies of the hg38 reference with random sequencing errors, drawn at
the per-base error rate implied by each base's own quality score. Read length,
insert size and base quality are sampled from distributions measured
in the real capture-panel BAMs, so the reads carry the same short,
adapter-trimmed, heavily overlapping fragment structure rather than an idealised
uniform layout:

| | Synthetic BAMs |
|---|---|
| Read length (median)  | 84 bp |
| Insert size (median) | 154 bp |
| MAPQ | 60 | 60 |
| Base quality (median) | 37 |
| Background non-reference rate | 2.8e-4 |
| Depth over target regions| ~18,100x |

At 18,000x it puts a handful of erroneous
reads under every position, so the caller is not being handed a noise-free
BAM. Depth is set above the real panel's so that the lowest designed
allele fractions (1e-4) are representable at all.

## Contents of `results/`

| | |
|---|---|
| `RESULTS.csv` | The scored table above |
| `calls/S<n>_<arm>.csv` | MNVista's passing calls, all 15 runs |
| `calls/S<n>_<arm>.falsepos.txt` | The false positives, listed |
| `logs/S<n>_<arm>.log` | MNVista run logs |
