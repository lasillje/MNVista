# MNVista

[![Build](https://github.com/lasillje/MNVista/actions/workflows/build_mnvista.yml/badge.svg)](https://github.com/lasillje/MNVista/actions/workflows/build.yml) [![License: MIT](https://img.shields.io/github/license/lasillje/MNVista)](LICENSE.md)

MNVista provides enhanced identification of multi-nucleotide variants from targeted short-read sequencing data. By utilising the detection of multi-nucleotide variants in cancers, the limit of detection (LOD) of minimal residual disease (MRD) can be increased.

MNVista was originally designed for use with B-cell lymphoma targeted short-read sequencing data.

## Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Building](#building)
- [Usage](#usage)
- [Recommended workflows](#recommended-workflows)
- [Output](#output)
- [License](#license)
- [Referencing](#referencing)

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


## Prerequisites

MNVista is written in C++23 and needs to be compiled from source. You will need:

- A C++23-capable compiler (e.g. GCC 13+, Clang 16+)
- [CMake](https://cmake.org/) 3.14 or newer
- [HTSlib](https://github.com/samtools/htslib) 1.19 (for sequence data parsing)
- OpenSSL development headers (`libssl-dev`)

**Installing the prerequisites (Ubuntu / Debian):**

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libhts-dev libssl-dev git
```

If you're on Windows, the simplest route is to install [WSL2](https://learn.microsoft.com/windows/wsl/install) with an Ubuntu distribution, then follow the Linux instructions above from within it.

## Building

Once the prerequisites above are installed, you can build MNVista as follows:

```bash
# 1. Clone the repository
git clone https://github.com/lasillje/MNVista.git
cd MNVista

# 2. Configure and build (this creates and uses a 'build' folder for you)
cmake -S . -B build
cmake --build build
```

If both commands complete without errors, the `MNVista` executable will be at `MNVista/bin/MNVista`. You can check that it worked by running:

```bash
./bin/MNVista --help
```

which should print the list of available options.

### Custom install locations

If HTSlib or OpenSSL are installed somewhere CMake can't find automatically (for example, a Conda environment or a manually built copy of HTSlib), you can point CMake directly at the headers and library files with:

```bash
cmake -S . -B build \
  -DHTSLIB_INCLUDE_DIR=/path/to/htslib/include \
  -DHTSLIB_LIBRARY=/path/to/htslib/lib/libhts.so \
  -DOPENSSL_ROOT_DIR=/path/to/openssl
cmake --build build
```

## Usage

MNVista expects at least 3 mandatory inputs:

- An indexed BAM file
- A VCF file sorted in ascending order
- The name of your desired output directory

The simplest run using all default parameters can then be done by the following:

```bash
MNVista IN_BAM IN_VCF OUT_DIR
```

This will output the results in OUT_DIR. By default, the name of the output files is 'results'.

A run using a few common options, e.g. naming the output, using 8 threads, and requiring a minimum read mapping quality of 30, would look like:

```bash
MNVista IN_BAM IN_VCF OUT_DIR -O my_sample -T 8 -Q 30
```

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
| -Q, --read-quality | Minimum read mapping quality |
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
| -K, --skip-filtered | Skip saving and outputting filtered MNVs |  
| -J, --min-jaccard | Minimum jaccard index for an MNV |
| -P, --bayes-prior-mnv | Bayes prior |

Note that any option regarding input SNV VAF or VRD requires this field to be present in the input VCF. Currently, accepted names for these fields are "AF" or "VAF" for the SNV VAFs, and "MRD" or "VRD" for the variant supporting read depth. In case these fields are not present, these options will be ignored.

For a comprehensive list of run options and default values, use ```MNVista -h```

## Recommended workflows

For the analysis of baseline samples (taken at time of diagnosis), we recommend to use a broad, unfiltered input VCF containing as many variants as possible (called at maximum sensitivity). Within MNVista, the minimum VRD (-S) of input SNVs can then be set to a value that best represents the VRD associated with the active disease. This heavily depends per disease and sequencing methods used. We used a value of 15 for coverages between ~15,000 - 22,000. Filter criteria should be quite stringent (B > 0.9, F > 0.5, J > 0.5) in order to call only high-confidence MNVs. Recommended maximum MNV size (-M) can be set to 2 (only SNV pairs) or higher (larger MNVs), however setting it to 0 (dynamic MNV size based on SNV count of current genomic window) or a large value ( > 10) might increase analysis time substantially. Lastly, if you are not interested in the filtered MNVs, the saving and outputting of these MNVs can be skipped (-K true), which also speeds up analysis time and decreases memory usage.

For follow-up samples (interim, end-of-treatment, etc.), the MNVista output VCF from the corresponding baseline should be used. This way, only SNVs associated with baseline MNVs will be considered, speeding up analysis time substantially. Filter criteria can be set very lax, as the presence of an MNV is now more important than the statistics. Minimum SNV VRD (-S) should be set to 1 for maximum sensitivity, and the maximum MNV size (-M) should be kept the same as in the baseline.

## Output

MNVista has 4 outputs:

- 'results.csv' containing MNVs that passed all criteria.
- 'results_filtered.csv' containing MNVs that failed one or more of the criteria.
- 'results.vcf' containing a list of all unique SNVs present in the results.csv, useful for keeping track of certain MNVs found in baseline samples across multiple timepoints.
- 'results.log' containing any log messages that were generated during a run, including the command used to generate these results.


| Output column | Description |
| --- | --- |
| WINDOW_ID | Internal ID of parent window |
| CHROM | Chromosome |
| MNV_NAME | Name of MNV, delimited by ':' and '-' |
| VAF | VAF of MNV |
| SD | Standard deviation of the MNV VAF |
| NUM_SUPPORTING | Read count of reads containing all MNV SNVs |
| NUM_COVERING | Read count of reads covering all MNV positions |
| PHI | Phi-coefficient of MNV |
| COEFFICIENT_VARIATION | CV of the MNV VAF |
| BAYESIAN | Bayesian probability of an MNV |
| JACCARD | Fraction of mutated reads containing all SNVs of MNV |
| DISCORDANT | Number of reads that contain only one SNV (per SNV) |
| INDIVIDUAL_MUTATED | VRD of each SNV |
| SIZE_MNV | Number of SNVs present in an MNV |
| DIST_MNV | Distance (bp) between first and last SNV of MNV |
| QUALITIES | Mean base quality of the reads supporting each SNV |
| LOG_ODDS | (Deprecated) Log-odds ratio for the MNV |
| DUPLICATE_FRACTION | Fraction of unique read position (start, end) combinations |
| MAX_SUPPORT_FRACTION | Fraction of reads stemming from the most common (start, end) combination |

## License
MNVista is licensed under the GPL-3.0 license. For more information, see LICENSE.md

## Referencing

_A citation for the MNVista paper will be added here once it is published._
