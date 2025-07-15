#pragma once

/*
    MNVista
    Copyright (C) 2025  Laurens Sillje

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
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