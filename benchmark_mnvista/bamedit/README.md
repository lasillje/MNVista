# Overview
BamEdit is a simple and fast tool for introducing spiked mutations in BAM files.
Given a .bed file and a .bam file, BamEdit will try to spike in the given mutations as accurately as possible.

# Recommended workflow
Create a .bed (tab/whitespace delimited) file with your desired mutations in the following format:

Chromosome  | Position | VAF | Haplotype 
------------- | ------------- | --- | ---
Chromosome (chr..) | Position of mutation | VAF of mutation | Haplotype ID


Haplotype ID is an integer which can be supplied in order to make sure certain mutations are spiked in the same reads.
For example, if there are 5 mutations on chr1 with a haplotype ID of '1' then BamEdit will spike these mutations into the same reads. This is naturally only supported for mutations that lie on the same chromosome. If no common reads between mutations of the same haplotype are found then the spiking cannot proceed and will be skipped.

Then, supply this bed file along with an indexed BAM file into BamEdit.

# Limitations
Make sure the reference fasta does not include any lower case (softmask) characters.
Currently, BamEdit is limited to only base substitutions. Additional support for indels will be added at a later point.

