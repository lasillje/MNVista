#!/usr/bin/env python3
"""
Generate clean synthetic paired-end BAMs spanning the 10 super-enhancer regions
used by the benchmark.

These BAMs are clean / contain no variants, every read is a copy of the hg38 reference
sequence apart from random sequencing errors, which are drawn at the per-base
error rate implied by that base's own quality score.

Read length, insert size and base quality are drawn from a distribution
measured from the BAMs used in the manuscript, so the synthetic reads carry largely the same structure.

The script is fully deterministic, each sample uses a fixed seed and rerunning this script
reproduces the BAMs exactly.

Usage:
    python make_synthetic_bams.py <hg38.fa> <regions_se.bed> <out_dir> [depth]

Writes <out_dir>/S1.bam ... S5.bam, coordinate-sorted and indexed.
"""

import os
import random
import sys
import pysam

READLEN_Q = [20, 26, 31, 35, 39, 42, 45, 48, 50, 53, 55, 58, 60, 62, 65, 67,
             70, 72, 74, 76, 78, 80, 82, 84, 87, 90, 92, 95, 98, 101, 105,
             108, 112, 116, 119, 123, 128, 133, 138, 142, 142, 142, 142, 142,
             142, 142, 142, 142, 142, 142]
INSERT_Q = [2, 27, 41, 51, 59, 67, 73, 79, 84, 89, 93, 97, 101, 105, 109, 113,
            117, 121, 125, 129, 132, 136, 139, 143, 146, 150, 154, 157, 160,
            164, 167, 172, 177, 185, 194, 205, 217, 230, 242, 253, 264, 275,
            286, 297, 309, 323, 342, 365, 395, 433]

BQ_MEAN, BQ_SD, BQ_MIN, BQ_MAX = 36.6, 1.9, 29, 40   # flat across the read
MAPQ = 60
MIN_READ = 25

SEEDS = {"S1": 20260901, "S2": 20260902, "S3": 20260903,
         "S4": 20260904, "S5": 20260905}

DEFAULT_DEPTH = 18000
FLANK = 500          # simulate fragments starting this far outside each region
_OTHER = {b: [x for x in "ACGT" if x != b] for b in "ACGT"}
_ERRP = [10 ** (-q / 10.0) for q in range(128)]


def read_regions(path):
    out = []
    for line in open(path):
        if line.strip():
            f = line.split()
            out.append((f[0], int(f[1]), int(f[2]),
                        f[3] if len(f) > 3 else "."))
    return out


def make_header(fa, regions):
    names = list(dict.fromkeys(r[0] for r in regions))
    names.sort(key=lambda c: int(c[3:]) if c[3:].isdigit() else 99)
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": fa.get_reference_length(c)} for c in names],
        "RG": [{"ID": "synthetic", "SM": "synthetic", "PL": "ILLUMINA"}],
    })


def expected_bases_per_fragment():
    tot = 0
    for r in READLEN_Q:
        for i in INSERT_Q:
            tot += 2 * min(r, max(i, MIN_READ))
    return tot / (len(READLEN_Q) * len(INSERT_Q))


def simulate_region(rng, fa, chrom, start, end, tid, depth, prefix):
    """Return [(pos, AlignedSegment)] for one region, sorted by position."""
    lo, hi = start - FLANK, end + FLANK
    ref = fa.fetch(chrom, lo, hi).upper()

    span = (hi - lo) - sum(INSERT_Q) / len(INSERT_Q)
    n_frag = int(round(depth * span / expected_bases_per_fragment()))

    choice, random_, gauss = rng.choice, rng.random, rng.gauss
    recs = []
    for i in range(n_frag):
        ins = max(choice(INSERT_Q), MIN_READ)
        rlen = min(choice(READLEN_Q), ins)
        if rlen < MIN_READ:
            continue
        fstart = rng.randint(lo, hi - ins - 1)
        qname = f"{prefix}:{i}"
        p1, p2 = fstart, fstart + ins - rlen

        for first, pos in ((True, p1), (False, p2)):
            seq = ref[pos - lo: pos - lo + rlen]
            if len(seq) < rlen or "N" in seq:
                continue
            quals = [max(BQ_MIN, min(BQ_MAX, int(round(gauss(BQ_MEAN, BQ_SD)))))
                     for _ in range(rlen)]
            # Substitute a base wherever its quality score says it is wrong.
            # At Q37 that is 2e-4 per base, which at 18,000x leaves a few errors
            out = None
            for j, q in enumerate(quals):
                if random_() < _ERRP[q]:
                    if out is None:
                        out = list(seq)
                    if out[j] in _OTHER:
                        out[j] = choice(_OTHER[out[j]])
            if out is not None:
                seq = "".join(out)

            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(33 + q) for q in quals))
            a.flag = 0
            a.is_paired = True
            a.is_proper_pair = True
            a.is_read1 = first
            a.is_read2 = not first
            a.is_reverse = not first
            a.mate_is_reverse = first
            a.reference_id = tid
            a.reference_start = pos
            a.next_reference_id = tid
            a.next_reference_start = p2 if first else p1
            a.template_length = ins if first else -ins
            a.mapping_quality = MAPQ
            a.cigarstring = f"{rlen}M"
            a.set_tag("RG", "synthetic")
            recs.append((pos, a))
    recs.sort(key=lambda x: x[0])
    return recs


def main():
    fa_path, bed_path, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    depth = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_DEPTH
    only = sys.argv[5] if len(sys.argv) > 5 else None
    os.makedirs(out_dir, exist_ok=True)

    fa = pysam.FastaFile(fa_path)
    regions = read_regions(bed_path)
    hdr = make_header(fa, regions)
    tids = {d["SN"]: i for i, d in enumerate(hdr.to_dict()["SQ"])}
    regions.sort(key=lambda r: (tids[r[0]], r[1]))

    for sample, seed in SEEDS.items():
        if only and sample != only:
            continue
        for ri, (chrom, start, end, gene) in enumerate(regions):
            out = os.path.join(out_dir, f"{sample}.part{ri:02d}.bam")
            if os.path.exists(out + ".done"):
                continue
            rng = random.Random(seed * 100 + ri)
            n = 0
            with pysam.AlignmentFile(out, "wb", header=hdr) as bam:
                for _, rec in simulate_region(rng, fa, chrom, start, end,
                                              tids[chrom], depth,
                                              f"{sample}:{gene}"):
                    bam.write(rec)
                    n += 1
            open(out + ".done", "w").write(str(n))
            print(f"{sample} {gene}: {n:,} reads", flush=True)


def merge_chunks(out_dir, sample, n_parts):
    """Concatenate the per-region chunks into one coordinate-sorted BAM."""
    parts = [os.path.join(out_dir, f"{sample}.part{i:02d}.bam")
             for i in range(n_parts)]
    out = os.path.join(out_dir, f"{sample}.bam")
    pysam.cat("-o", out, *parts)      # chunks are already in genomic order
    pysam.index(out)
    for p in parts:
        os.remove(p)
        os.remove(p + ".done")
    return out


if __name__ == "__main__":
    main()
