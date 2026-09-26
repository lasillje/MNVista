import sys
import os
import random
import statistics

VAFS = [0.001, 0.005, 0.01, 0.05, 0.1]
SPIKE_VAFS = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.1]

REGIONS_FILE = sys.argv[1]
HAPLOTYPE_SIZE = int(sys.argv[2])
OUTPUT_NAME = sys.argv[3]
NOISE_FRAC = float(sys.argv[4])
MAX_CLUSTER = int(sys.argv[5])

SNP_LINES = []
BED_LINES = []
TN_LINES = []
MNV_LINES = []

VCF_LINES = []

OCCUPIED_POSITIONS = set()

REGIONS = []

VAF_COUNTS = [0] * len(VAFS)
SPIKE_VAF_COUNTS = [0] * len(SPIKE_VAFS)

HAPLO_COUNTER = 0

SNV_COUNTER = 0
NOISE_COUNTER = 0
MNV_COUNTER = 0

random.seed()


def add_noise(noise_level, lines: list):

    NOISE_LINES = []
    
    count = 0
    for line in lines:
        s = line.strip().split()
        NOISE_LINES.append(f"{s[0]} {s[1]} {float(s[2]) + (float(s[2]) * random.uniform(-noise_level, noise_level))} {s[3]}")
        count+=1
    write_lines(NOISE_LINES, OUTPUT_NAME, f"_NOISE_{noise_level}", "bed", False)

def gen_single_snp(region: list, pos: int):
    if(pos not in OCCUPIED_POSITIONS):
        SNP_LINES.append(f"{region[0]} {pos} {pos}  0.5")
        OCCUPIED_POSITIONS.add(pos)


def gen_cluster(region: list, start_pos, end_pos, sz):
    global HAPLO_COUNTER
    global MNV_COUNTER
    global SNV_COUNTER
    global NOISE_COUNTER

    vaf = random.choice(VAFS)
    if(sz == 1):
        vaf = random.choice(SPIKE_VAFS)
    size = sz
    positions = random.sample(range(start_pos, end_pos), size)
    positions.sort()

    if sz > 1:
        if positions[0] == positions[1]:
            positions[1] += 1

    gen = True
    for i in range(0, size):
        if(positions[i] in OCCUPIED_POSITIONS):
            gen = False
            break
    if gen:
        for i in range(0, size):
            BED_LINES.append(f"{region[0]} {positions[i]} {vaf} {HAPLO_COUNTER}")
            OCCUPIED_POSITIONS.add(positions[i])
        HAPLO_COUNTER += 1
        if sz == 1:
            NOISE_COUNTER += 1
            SNV_COUNTER += 1
            SPIKE_VAF_COUNTS[SPIKE_VAFS.index(vaf)] += 1
        else:
            MNV_COUNTER += 1
            SNV_COUNTER += 2
            VAF_COUNTS[VAFS.index(vaf)] += size

    if sz > 1 and gen:
        mnv_name = "-".join(str(x) for x in positions)
        MNV_LINES.append(f"{region[0]}:{mnv_name}")

def gen_clusters_in_region(region: list, size):
    if size > 1:
        cluster_dif = region[2] - region[1]
        max_clusters = MAX_CLUSTER
        if(max_clusters > 0):
            for i in range(0, int(max_clusters)):
                start_pos = random.randint(region[1], region[2] - HAPLOTYPE_SIZE - 1)
                if(start_pos + HAPLOTYPE_SIZE < region[2]):
                    gen_cluster(region, start_pos, start_pos + HAPLOTYPE_SIZE, size)
    else:
        cluster_dif = region[2] - region[1]
        cluster_dif = cluster_dif * float(NOISE_FRAC)
        if(int(cluster_dif) == 0):
            return
        if(cluster_dif <= 1):
            cluster_dif = 3
        for i in range(0, int(cluster_dif)):
            pos = random.randint(region[1], region[2] - 1)
            gen_cluster(region, pos, pos+1, 1)

def gen_true_positives(bed_regions: list):
    for region in bed_regions:
        gen_clusters_in_region(region, 2)
        gen_clusters_in_region(region, 1)

    write_lines(BED_LINES, OUTPUT_NAME, "", "bed", True)
    write_lines(MNV_LINES, OUTPUT_NAME, "mnv", "csv", False)
    
    add_noise(0.1, BED_LINES)
    
    OCCUPIED_POSITIONS.clear()

def gen_clusters_bed(bed_regions: list):
    gen_true_positives(bed_regions)

def sort_func(x):
    splitted_lines = x.split(' ')
    chrom_num = splitted_lines[0].lower().lstrip('chr')
    num = 0
    if chrom_num == 'X':
        num = 9998
    elif chrom_num == 'Y':
        num = 9999
    else:
        num = int(chrom_num)
    return (num, int(splitted_lines[1]))

def write_lines(lines: list, outname: str, suffix: str, format: str, do_sort):
    if do_sort:
        lines.sort(key = sort_func)
    out_name = f"{outname}_{suffix}.{format}"
    if not suffix:
        out_name = f"{outname}.{format}"
    with open(out_name, "w") as f:
        for line in lines:
            f.write(line + "\n")
    with open(f"{outname}.stats", "w") as f2:
        mean_counts = statistics.mean(VAF_COUNTS)
        sd_counts = statistics.stdev(VAF_COUNTS)
        noise_mean = statistics.mean(SPIKE_VAF_COUNTS)
        noise_sd = statistics.stdev(SPIKE_VAF_COUNTS)
        f2.write("total;mean;sd;noise;noise_0;noise_1;noise_2;noise_3;noise_4;noise_5;noise_mean;noise_sd\n")
        f2.write(f"{SNV_COUNTER-NOISE_COUNTER};{mean_counts};{sd_counts};{NOISE_COUNTER};{SPIKE_VAF_COUNTS[0]};{SPIKE_VAF_COUNTS[1]};{SPIKE_VAF_COUNTS[2]};{SPIKE_VAF_COUNTS[3]};{SPIKE_VAF_COUNTS[4]};{SPIKE_VAF_COUNTS[5]};{noise_mean};{noise_sd}\n")

    

with open(REGIONS_FILE, "r") as f:
    for line in f:
        if line.startswith("chr"):
            s = line.strip().split()
            REGIONS.append([s[0], int(s[1]), int(s[2])])
    print(f"Found {len(REGIONS)} regions in bed file.")

gen_clusters_bed(REGIONS)

print(f"Generated {SNV_COUNTER} SNVs.")
print(f"Generated {MNV_COUNTER} MNVss.")
print(f"Generated {NOISE_COUNTER} noise SNVs.")
print("VAF counts: ")
print(f"    0.001   {VAF_COUNTS[0]}")
print(f"    0.005   {VAF_COUNTS[1]}")
print(f"    0.01    {VAF_COUNTS[2]}")
print(f"    0.05    {VAF_COUNTS[3]}")
print(f"    0.1     {VAF_COUNTS[4]}")
