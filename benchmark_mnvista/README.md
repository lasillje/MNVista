# benchmark_mnvista

Everything needed to reproduce the synthetic spike-in benchmark reported in the
manuscript: the spiking tool, the spike-in generator, the target
regions, and the complete ground-truth call set for all five samples and all
three noise arms.

The sequencing data itself (the unspiked patient BAMs) is not redistributed
here. Given a BAM, the files in `truth_sets/` are sufficient to rebuild the
spiked BAMs, and to re-score any caller against the same truth set.

```
benchmark_mnvista/
├── bamedit/                 the spiking tool (C++, MIT)
├── spike_in_design/         generator + target regions
├── truth_sets/              the spiked variants (BED / VCF / MNV intervals)
└── validation/              a complete rerun of the benchmark, with results
```

`validation/` contains a self-contained end-to-end rerun of the synthetic benchmark as in the manuscript.
Instead of patient BAMs, clean synthetic BAMs are generated over the target regions, the deposited truth sets are spiked into
them, MNVista is called with the manuscript's settings, and the calls are scored
against the truth. All 8,008 spiked MNVs were recovered, with no misses and 41
false positives at a minimum supporting read of 1, dropping to a single false positive at a minimum supporting read of 3 (recall 1.000, precision 0.995). See `validation/README.md`.

---

## 1. `bamedit/` (tool for spiking MNVs)

`BamEdit` is the tool referred to in the Methods as the modified BAMSurgeon-style
spike-in workflow. It takes an indexed BAM and a BED of `chrom / pos / VAF / haplotype-ID` and writes the SNVs into the reads. Variants sharing a
haplotype ID are spiked into the same reads, to create an actual MNV. Only base substitutions are supported.

Build:

```bash
cd bamedit && mkdir build && cd build
cmake .. && make
```

Run:

```bash
bamedit <spike_ins.bed> <reference.fa> <input.bam> <output_dir> -O <name> -S <seed>
samtools index <output_dir>/<name>.bam
```

This writes `<output_dir>/<name>.bam` (the spiked BAM) and `<output_dir>/<name>.vcf`
(the variants as actually written, with realised `AF`, `VRD`, `DP` and `HAPLO`).
This VCF is the candidate list MNVista is then run with. `-S/--seed` makes read
selection deterministic, the runs reported here used the default, `-S 1`.

BamEdit currently takes the reference base literally, so at a lowercase (softmasked) position it can pick the uppercase
form of the same base as the alternate allele and write no mutation at all.
So, first uppercase the FASTA before use (if not done already):

```bash
awk '/^>/ {print; next} {print toupper($0)}' hg38.fa > hg38.upper.fa
samtools faidx hg38.upper.fa
```

See `bamedit/README.md` for the full interface.

---

## 2. `spike_in_design/` (spike in generation)

| File | Usecase |
|---|---|
| `bedgen.py` | The generator that produced every BED in `truth_sets/beds/` |
| `regions_se.bed` | The 10 super-enhancer loci the spike-ins were placed in (same regions used in benchmark) |
| `se_regions.csv` | The same 10 loci with gene names (Supplementary Table S1) |

```bash
python bedgen.py regions_se.bed 100 <output_prefix> <noise_fraction> <mnvs_per_region>
```

- argument 2 (`100`) is the haplotype window: the SNVs of one MNV are drawn
  from a 100 bp window, matching the manuscript method.
- `noise_fraction` sets how densely single-SNV noise is seeded across each
  region (differs per the 0%, 1% and 10% noise analysis)
- Every MNV in this benchmark consists of exactly **2 SNVs**.
- The script also writes a `__NOISE_0.1` copy of the BED in which each VAF is
  jittered by up to ±10%, and an `_mnv.csv` listing the MNV intervals.

### Note on random seeds

The used `bedgen.py` calls `random.seed()` with no argument, so the designs were drawn
from an unrecorded seed and the exact seed cannot be recovered.
This does not affect reproducibility of the benchmark, as the generator's output
is the ground truth, and the complete output is deposited here in
`truth_sets/`. Re-running `bedgen.py` produces a new, equally valid design.

`validation/` demonstrates this, the whole benchmark is rerun from the deposited files on freshly generated synthetic BAMs, and recovers
every one of the 8,008 spiked MNVs.

---

## 3. `truth_sets/` (the ground truth)


| Directory | Contents |
|---|---|
| `beds/` | `S<n>_<arm>.bed` exactly the file passed to BamEdit: `chrom  pos  VAF  haplotype_id` |
| `vcfs/` | `S<n>_<arm>.vcf` the same variants as VCF, with the realized `AF`, variant read depth `VRD`, total depth `DP` and `HAPLO` (haplotype ID) measured in the spiked BAM |
| `mnv/` | `S<n>_<arm>_mnv.csv` the MNV-level truth set, one `chrom:pos1-pos2` interval per spiked MNV |
| `manifest.csv` | Per-file counts, plus the original filenames each file was derived from |


**Arms.** `noise0` = no added noise (MNVs only); `noise1` and `noise10` add
single-SNV background noise at increasing density. The three arms are
independent draws, not the same MNVs with noise layered on, so the MNV sets
differ between arms. For `noise1` and `noise10` the BED that was actually spiked
is the ±10% VAF-jittered one, which is why those VAFs are not round numbers.
`noise0` was spiked at the exact designed VAFs.

**Counts.**

| Sample | Arm | Spiked SNVs | MNVs | Noise SNVs |
|---|---|---|---|---|
| S1 | noise0 | 1,060 | 530 | 0 |
| S1 | noise1 | 1,178 | 544 | 90 |
| S1 | noise10 | 1,917 | 531 | 855 |
| S2 | noise0 | 1,060 | 530 | 0 |
| S2 | noise1 | 1,163 | 537 | 89 |
| S2 | noise10 | 1,897 | 523 | 851 |
| S3 | noise0 | 1,066 | 533 | 0 |
| S3 | noise1 | 1,158 | 535 | 88 |
| S3 | noise10 | 1,943 | 546 | 851 |
| S4 | noise0 | 1,090 | 545 | 0 |
| S4 | noise1 | 1,164 | 536 | 92 |
| S4 | noise10 | 1,888 | 521 | 846 |
| S5 | noise0 | 1,054 | 527 | 0 |
| S5 | noise1 | 1,149 | 532 | 85 |
| S5 | noise10 | 1,924 | 538 | 848 |
| **Total** | | **20,711** | **8,008** | **4,695** |

---

## 4. Reproducing the benchmark

See validation/README.md for a comprehensive guide on re-running the benchmark.

