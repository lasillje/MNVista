
# MNVista

MNVista provides enhanced identification of multi-nucleotide variants from targeted short-read sequencing data. By utilising the detection of multi-nucleotide variants in cancers, the limit of detection (LOD) of minimal residual disease (MRD) can be increased.

MNVista was originally designed for use with B-cell lymphoma targeted short-read sequencing data

## Overview

The following chapters will make use of various abbreviations, if you are not familiar with them please refer to the following table.

| Option | Description |
| --- | --- |
| MNV | Multi-nucleotide variant |
| SNV | Single nucleotide variant |
| VAF | Variant allele frequency |
| VRD | Variant read depth |
| CV | Coefficient of variation |
| bp | Base Pair |


## Dependencies

MNVista is written in C++23 using HTSlib 1.19 for sequence data parsing and argparse for commandline functionality.
  - HTSlib
  - OpenSSL (libssl-dev)

## Building

To build MNVista, run the following commands in order:

```bash
# Clone or download the MNVista repo:
git clone https://github.com/lasillje/MNVista.git
cd MNVista

# Create build directory
mkdir build

# OPTIONAL: Configure dependency paths in CMake
cmake .. \
  -DHTSLIB_ROOT=/path/to/htslib \
  -DOpenSSL_ROOT_DIR=/path/to/openssl

# Build
cmake -B build && cmake --build build
```
The final executable can be found in MNVista/bin/

## Usage

MNVista expects atleast 3 mandatory inputs:

- An indexed BAM file
- A VCF file sorted in ascending order
- The name of your desired output directory

The simplest run using all default parameters can then be done by the following:

```bash
MNVista IN_BAM IN_VCF OUT_DIR
```

This will output the results in OUT_DIR. By default, the name of the output files is 'results'.

An overview of currently available run options with a short description:

| Option | Description |
| --- | --- |
| input_bam | Path to input BAM file|
| input_vcf | Path to input VCF file |
| output_dir | Path to output directory |
| -h, --help | List all available options |
| -v, --version | Show MNVista version |
| -V, --verbose | Enable verbose logging |
| -O, --output-name | Sets output file name | 
| -M, --max-mnv-size | Sets max SNVs in an MNV |
| -Q, --read-quality | Minimum read quality |
| -T, --threads | Amount of threads to be used|
| -Y, --min-vaf-snv | Minimum VAF of input SNVs |
| -A, --min-vaf-mnv | Minimum VAF of output MNVs |
| -S, --min-vrd-snv | Minimum VRD of input SNVs |
| -N, --min-vrd-mnv | Minimum VRD of output MNVs |
| -Z, --max-vaf-snv | Maximum VAF of input SNVs |
| -B, --min-bayesian | Minimum Bayes probability for MNVs|
| -F, --min-phi | Minimum Phi-coefficient for MNVs|
| -C, --black-list | List of MNVs to be ignored|
| -R, --read-length | Expected length of reads in data (bp)|  
| -E, --bayes-error-freq | Minimum expected reliable SNV VAF |
| -H, --bayes-mnv-freq | Minimum expected MNV VAF |
| -P, --bayes-prior-mnv | Bayes prior for an MNV being real |


For a comprehensive list of run options and default values, use ```MNVista -h```

## Output

MNVista has 4 outputs:

- 'results.csv' containing MNVs that have a higher Bayesian probability and Phi-coefficient than the minimum specified value.
- 'results_filtered.csv' containing MNVs that failed the previous criteria.
- 'results.vcf' containing a list of all unique SNVs present in the results.csv, useful for keeping track of certain MNVs found in baseline samples across multiple timepoints.
- 'results.log' containing any log messages that were generated during a run.


| Output column | Description |
| --- | --- |
| WINDOW_ID | Internal ID of parent window |
| CHROM | Chromosome |
| MNV_NAME | Name of MNV, delimited by ':' and '-'
| VAF | VAF of MNV |
| SD | Standard deviation of the MNV VAF |
| NUM_SUPPORTING | Read count of reads containing all MNV SNVs |
| NUM_COVERING | Read count of reads covering all MNV positions |
| PHI | Phi-coefficient of MNV |
| COEFFICIENT_VARIATION | CV of the MNV VAF |
| BAYESIAN | Bayesian probability of an MNV |
| JACCARD | Fraction of mutated reads containing all SNVs of MNV|
| DISCORDANT | Number of reads that contain only one SNV (per SNV)| 
| INDIVIDUAL_MUTATED| VRD of each SNV|
|SIZE_MNV | Number of SNVs present in an MNV|
| DIST_MNV | Distance (bp) between first and last SNV of MNV |


## Referencing
