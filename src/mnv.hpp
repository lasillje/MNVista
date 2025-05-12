#pragma once

#define mnv_window std::vector<mnv>
#define snv_window std::vector<snv*>

#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <fstream>

enum MNV_RESULT : unsigned char
{
    MNV_SUCCESS = 0,
    MNV_NO_SHARED_READS = 1,
    MNV_FAILED_FILTERS = 2
};

struct run_params
{
    std::string in_bam;
    std::string in_vcf;
    std::string out_path;
    std::string out_name;
    std::string blacklist_path;
    std::ofstream log_stream;
    double odds_ratio;
    double min_bayesian;
    double min_phi;
    int num_threads;
    int min_read_quality;
    int window_size;
    int mnv_size;
    int min_mrd_snv;
    int min_mrd_mnv;
    int read_length;
    float min_snv_vaf;
    float max_snv_vaf;
    float min_vaf;
    bool verbose;
    bool zero_indexed;
};

struct snv
{
    std::vector<unsigned int> supporting_hashes;
    std::vector<unsigned int> covering_hashes;
    std::vector<unsigned int> base_qualities;
    std::string chrom_name;
    double phred_qual;
    int chrom_id;
    unsigned int pos;
    unsigned int mrd;
    unsigned int dp;
    float vaf;
    char ref;
    char alt;
    char loaded_reads;
};

struct mnv
{
    std::string name;
    std::vector<unsigned short> discordant;
    snv_window variants;
    double odds_ratio;
    double odds_phi;
    double bayesian_prob;
    float vaf;
    float frac;
    float sd;
    float rsd;
    unsigned int window_id;
    unsigned int num_sup;
    unsigned int num_cov;
    unsigned char size;
    unsigned char result;
};

struct mnv_hash
{
  size_t operator()(const mnv& m) const
  {
    return std::hash<std::string>()(m.name);
  }
};

extern run_params settings;