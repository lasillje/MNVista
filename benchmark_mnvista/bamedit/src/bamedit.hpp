#include <string>

struct snv
{
    std::string chrom;
    int pos;
    char ref, alt;
    double vaf;
    int vrd;
    int dp;
    int haplotype;
    int tid;
};

