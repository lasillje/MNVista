#!/usr/bin/env python3
"""
Score MNVista calls against the deposited MNV truth sets.

A call is a true positive when the set of genomic positions it is built from is
exactly a spiked haplotype.

Because the candidate VCF handed to MNVista is BamEdit's own output, every
candidate SNV is a spiked SNV. A false positive is therefore a purely pairing error:
MNVista joining two SNVs that were never (explicitly) spiked onto the same reads.
Ofcourse, by chance, SNVs can still be spiked into the same read, which would be the main source of false positives here.

Usage:
    python score_calls.py <truth_dir> <calls_dir> <out.csv>

truth_dir  ../truth_sets/mnv          (S<n>_<arm>_mnv.csv)
calls_dir  directory of MNVista .csv  (S<n>_<arm>.csv)
"""

import csv
import glob
import os
import sys

def truth_keys(path):
    """{(chrom, (pos, ...))} from a truth MNV csv (chrom:pos1-pos2)."""
    out = set()
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        chrom, rest = line.split(":", 1)
        out.add((chrom, tuple(sorted(int(x) for x in rest.split("-")))))
    return out


def call_keys(path):
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter=";"):
            name = row["MNV_NAME"]
            chrom, coords = name.split(":")[0], name.split(":")[1]
            key = (chrom, tuple(sorted(int(x) for x in coords.split("-"))))
            out[key] = row
    return out


def main():
    truth_dir, calls_dir, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
    rows = []
    for cp in sorted(glob.glob(os.path.join(calls_dir, "S*.csv"))):
        name = os.path.basename(cp)[:-4]
        if name.endswith("_filtered"):
            continue
        tp_path = os.path.join(truth_dir, f"{name}_mnv.csv")
        if not os.path.exists(tp_path):
            continue
        truth = truth_keys(tp_path)
        calls = call_keys(cp)
        tp = truth & set(calls)
        fp = set(calls) - truth
        fn = truth - set(calls)
        prec = len(tp) / len(calls) if calls else float("nan")
        rec = len(tp) / len(truth) if truth else float("nan")
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        sample, arm = name.split("_", 1)
        rows.append(dict(sample=sample, arm=arm, truth_mnvs=len(truth),
                         calls=len(calls), tp=len(tp), fp=len(fp), fn=len(fn),
                         precision=round(prec, 6), recall=round(rec, 6),
                         f1=round(f1, 6)))
        # keep the misses for inspection, most likely does nothing
        with open(os.path.join(calls_dir, f"{name}.missed.txt"), "w") as fh:
            for c, p in sorted(fn):
                fh.write(f"{c}:{'-'.join(map(str, p))}\n")
        # also keep falsepositives overview
        with open(os.path.join(calls_dir, f"{name}.falsepos.txt"), "w") as fh:
            for c, p in sorted(fp):
                fh.write(f"{c}:{'-'.join(map(str, p))}\n")

    tot = {k: sum(r[k] for r in rows) for k in
           ("truth_mnvs", "calls", "tp", "fp", "fn")}
    p = tot["tp"] / tot["calls"] if tot["calls"] else 0
    r = tot["tp"] / tot["truth_mnvs"] if tot["truth_mnvs"] else 0
    rows.append(dict(sample="ALL", arm="all", **tot,
                     precision=round(p, 6), recall=round(r, 6),
                     f1=round(2 * p * r / (p + r), 6) if p + r else 0.0))

    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for row in rows:
        print("  ".join(f"{k}={v}" for k, v in row.items()))


if __name__ == "__main__":
    main()
