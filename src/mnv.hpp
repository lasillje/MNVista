#pragma once

/*
MIT License

Copyright (c) 2025 Laurens Sillje

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


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
    double bayes_freq;
    double bayes_haplo;
    double bayes_prior;
    double jaccard;
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
    bool skip_filtered;
};

struct read
{
    std::string read_name;
    unsigned char quality;

    bool operator < (read const& other) const noexcept
    {
      return read_name < other.read_name;
    }

    bool operator == (read const& other) const noexcept
    {
      return read_name == other.read_name;
    }
};

struct snv
{
    std::vector<read> supporting_hashes;
    std::vector<read> covering_hashes;
    //std::vector<unsigned int> base_qualities;
    std::string chrom_name;
    double base_qual_sum;
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
    std::vector<float> qualities;
    std::vector<float> discordant_qualities;
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