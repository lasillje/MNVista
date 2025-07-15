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

#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <sstream>
#include <filesystem>
#include <chrono>
#include <functional>
#include <mutex>

#include "mnv.hpp"
#include "stats.hpp"
#include "window.hpp"
#include "utils.hpp"

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/thread_pool.h"
#include "htslib/kfunc.h"

#include "argparse/argparse.hpp"

bool operator==(const mnv& lhs, const mnv& rhs)
{
    return lhs.name == rhs.name;
};

std::unordered_set<std::string> pair_blacklist;

std::set<int> chrom_contigs;
std::set<std::string> vcf_contigs;

std::vector<snv> all_snvs;
std::vector<snv_window> all_snv_windows;

mnv_window* all_mnv_windows;
mnv_window* all_mnv_filtered;

int total_windows = 0;
bool blacklist_enable = false;

run_params settings;

/*
    Reads input VCF list and stores it into variant_list
*/
void read_vcf(std::vector<snv> &variant_list, std::string in_vcf, std::string in_bam)
{

    samFile* bam = sam_open(in_bam.c_str(), "r");
    if(bam == nullptr)
    {
      return;
    }
    
    bam_hdr_t* bam_hdr = sam_hdr_read(bam);
    if(bam_hdr == nullptr)
    {
      return;
    }
    
    hts_idx_t* bam_idx = sam_index_load(bam, in_bam.c_str());
    if(!bam_idx)
    {
      return;
    }


    htsFile *vcf = bcf_open(in_vcf.c_str(), "r");
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);
    bcf1_t *rec = bcf_init();

    while (bcf_read(vcf, hdr, rec) == 0)
    {
        bcf_unpack(rec, BCF_UN_STR);

        snv v{};
        v.loaded_reads = 0;
        v.pos = rec->pos;

        v.chrom_name = std::string(bcf_hdr_id2name(hdr, rec->rid));
        int tid = bam_name2id(bam_hdr, v.chrom_name.c_str());

        if(tid < 0)
        {
            continue;
        }

        if (std::strlen(rec->d.allele[0]) != 1 || std::strlen(rec->d.allele[1]) != 1)
        {
            if (settings.verbose)
            {
                log_info("Skipped indel at chr" + v.chrom_name + ":" + std::to_string(v.pos) + "\n");
            }
            continue;
        }

        v.ref = rec->d.allele[0][0];
        v.alt = rec->d.allele[1][0];

        float *vafs = nullptr;
        int vaf_count = 0;

        bool found_vaf_field = true;
        if ((bcf_get_info_float(hdr, rec, "AF", &vafs, &vaf_count) >= 0 || bcf_get_info_float(hdr, rec, "VAF", &vafs, &vaf_count) >= 0 ) && vaf_count > 0)
        {
            v.vaf = vafs[0];
            free(vafs);
        }
        else
        {
            found_vaf_field = false;
            v.vaf = 1.0f;
        }

        int* mrds = nullptr;
        int mrd_count = 0;
        bool found_mrd_field = true;
        if( (bcf_get_info_int32(hdr, rec, "MRD", &mrds, &mrd_count) >= 0 || bcf_get_info_int32(hdr, rec, "VRD", &mrds, &mrd_count) >= 0  ) && mrd_count > 0)
        {
            v.mrd = mrds[0];
            free(mrds);
        } else
        {
            v.mrd = 1;
            found_mrd_field = false;
        }

        if (found_vaf_field && (v.vaf < 0.0f || v.vaf < settings.min_snv_vaf || v.vaf > settings.max_snv_vaf))
        {
            if (settings.verbose)
            {
                log_info("Skipped variant at " + v.chrom_name + ":" + std::to_string(v.pos) + " due to VAF of " + std::to_string(v.vaf) + "\n");
            }
            continue;
        }

        if (found_mrd_field && v.mrd < settings.min_mrd_snv)
        {
            if (settings.verbose)
            {
                log_info("Skipped variant at " + v.chrom_name + ":" + std::to_string(v.pos) + " due to VRD of " + std::to_string(v.mrd) + "\n");
            }
            continue;
        }

        vcf_contigs.insert(v.chrom_name);
        chrom_contigs.insert(tid);
        v.chrom_id = tid;

        variant_list.push_back(v);
    }

    log_info("Found " + std::to_string(variant_list.size()) + " total variants\n");

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);

    sam_hdr_destroy(bam_hdr);
    hts_idx_destroy(bam_idx);
    sam_close(bam);
}

/*
    Tests any number of SNVs whether they make up an MNV or not
    Statistical tests are only performed for doublets
*/
MNV_RESULT test_snv(const std::vector<snv *> &s, int num_variants, mnv *out_mnv)
{
    out_mnv->size = (unsigned char)num_variants;
    out_mnv->sd = 0;

    std::sort(s[0]->covering_hashes.begin(), s[0]->covering_hashes.end());
    std::sort(s[0]->supporting_hashes.begin(), s[0]->supporting_hashes.end());

    std::set<read> intersect_sup(s[0]->supporting_hashes.begin(), s[0]->supporting_hashes.end());
    std::set<read> intersect_cov(s[0]->covering_hashes.begin(), s[0]->covering_hashes.end());

    for (int i = 1; i < num_variants; i++)
    {
        std::sort(s[i]->covering_hashes.begin(), s[i]->covering_hashes.end());
        std::sort(s[i]->supporting_hashes.begin(), s[i]->supporting_hashes.end());

        std::set<read> temp_sup;
        std::set<read> temp_cov;

        std::set_intersection(intersect_sup.begin(), intersect_sup.end(), s[i]->supporting_hashes.begin(), s[i]->supporting_hashes.end(), std::inserter(temp_sup, temp_sup.begin()));
        std::set_intersection(intersect_cov.begin(), intersect_cov.end(), s[i]->covering_hashes.begin(), s[i]->covering_hashes.end(), std::inserter(temp_cov, temp_cov.begin()));

        intersect_sup = temp_sup;
        intersect_cov = temp_cov;
    }

    if (num_variants == 2 && intersect_sup.empty())
    {
        return MNV_NO_SHARED_READS;
    }

    int numMutated = intersect_sup.size();
    int numCov = intersect_cov.size();

    int numTotalSolo = 0;

    for (int i = 0; i < num_variants; i++)
    {
        std::set<read> onlyMyReads(intersect_cov.begin(), intersect_cov.end());
        std::set<read> myReads(s[i]->supporting_hashes.begin(), s[i]->supporting_hashes.end());

        std::erase_if(onlyMyReads, [&myReads](read r)
                      { return !myReads.contains(r); });

        for (int j = 0; j < num_variants; j++)
        {
            if (j == i)
                continue;
            std::set<read> otherReads(s[j]->supporting_hashes.begin(), s[j]->supporting_hashes.end());
            std::erase_if(onlyMyReads, [&otherReads](read r)
                          { return otherReads.contains(r); });
        }

        float sum_quality = 0.0;
        for(auto& read : onlyMyReads)
        {
            int q = (int) read.quality;
            float base_qual = (float)q / 10.0f;
            sum_quality += base_qual;
        }

        float concordant_quality = 0.0f;
        for(auto& read1 : intersect_sup)
        {
            for(auto& read2 : s[i]->supporting_hashes)
            {
                if(read1.read_name == read2.read_name)
                {
                    int q = (int)read2.quality;
                    concordant_quality += (float)q / 10.0f;
                }
            }
        }

        out_mnv->qualities.push_back(concordant_quality);
        out_mnv->discordant_qualities.push_back(sum_quality);
        out_mnv->discordant.push_back(onlyMyReads.size());
        numTotalSolo += onlyMyReads.size();
    }

    int mutTotal = numTotalSolo + numMutated;

    float totalFrac = mutTotal > 0 ? ((float)numMutated / (float)mutTotal) : 0.0f;
    out_mnv->vaf = numCov > 0 ? ((float)numMutated / (float)numCov) : 0.0f;
    
    out_mnv->num_sup = numMutated;
    out_mnv->num_cov = numCov;;

    float mean = 0.0f;
    out_mnv->sd = vaf_sd(s, num_variants, &mean);
    out_mnv->rsd = out_mnv->sd / mean;


    if(num_variants == 2)
    {
        out_mnv->frac = totalFrac;
        if(out_mnv->discordant.size() > 1)
        {
            int none = out_mnv->num_cov - out_mnv->num_sup - out_mnv->discordant[0] - out_mnv->discordant[1];
            out_mnv->odds_ratio = test_odds(out_mnv->num_sup, out_mnv->discordant[0], out_mnv->discordant[1], none);
            out_mnv->odds_phi = test_phi(out_mnv->num_sup, out_mnv->discordant[0], out_mnv->discordant[1], none);
            out_mnv->bayesian_prob = test_bayesian(s[0], s[1], out_mnv, out_mnv->num_sup, out_mnv->discordant[0], out_mnv->discordant[1], none, settings.bayes_freq, settings.bayes_haplo, settings.bayes_prior);
            if(settings.verbose)
            {
                std::stringstream debugstr;
                debugstr << out_mnv->bayesian_prob << "\n";
                std::cerr << debugstr.str();
            }
        }

        if (out_mnv->frac < settings.jaccard || out_mnv->odds_phi < settings.min_phi || out_mnv->odds_ratio < settings.odds_ratio || out_mnv->bayesian_prob < settings.min_bayesian || out_mnv->num_sup < settings.min_mrd_mnv || out_mnv->vaf < settings.min_vaf)
        {
            return MNV_FAILED_FILTERS;
        }
    }

    return MNV_SUCCESS;
}

/*
    Checks whether an MNV already contains a SNV based on chromosome and position
*/
bool mnv_contains_snv(mnv *m, snv *s)
{
    for (int i = 0; i < m->variants.size(); i++)
    {
        if (s->chrom_id == m->variants[i]->chrom_id && s->pos == m->variants[i]->pos)
        {
            return true;
        }
    }
    return false;
}

/*
    Tests a doublet for MNV creation
    Cached pairs will use an early out to skip redundant testing
*/
int test_doublet(snv *a, snv *b, mnv *out_mnv, std::unordered_set<mnv, mnv_hash> pair_cache)
{
    if(a->pos == b->pos || std::abs((int)a->pos - (int)b->pos) > settings.read_length) return 2;

    mnv temp = {};
    temp.variants.reserve(2);
    temp.variants.push_back(a);
    temp.variants.push_back(b);

    temp.name = name_mnv(temp.variants, true);

    if(blacklist_enable)
    {
        if(pair_blacklist.find(temp.name) != pair_blacklist.end())
        {
            std::stringstream bstr;
            bstr << "Found blacklisted pair at " << temp.name << ", skipping.";
            log_info(bstr.str());
            return 1;
        }
    }

    auto cached_pair = pair_cache.find(temp);

    if (cached_pair != pair_cache.end())
    {
        *out_mnv = *cached_pair;
        return out_mnv->result;
    }

    int res = test_snv(temp.variants, 2, &temp);
    temp.result = res;
    
    *out_mnv = temp;

    return res;
}

/*
    Makes a doublet from two SNVs
    out_result pointer is used to store the result of the MNV creation (failed/success)
    If the doublet is already cached, then the cached result is used instead
*/
mnv make_doublet(snv *a, snv *b, int *out_result, std::unordered_set<mnv, mnv_hash> &pair_cache)
{
    mnv temp = {};
    *out_result = test_doublet(a, b, &temp, pair_cache);
    return temp;
}

/*
    Makes the next phase size of an MNV
    For example, an MNV of size two turns into an MNV of size 3, if the additional SNV can make pairs with every already present SNV
*/
mnv make_next_phase(mnv *m, snv *extra_snv, int *out_result, std::unordered_set<mnv, mnv_hash> &pair_cache)
{
    bool can = true;
    for (int i = 0; i < m->variants.size(); i++)
    {
        mnv temp = {};
        int res = test_doublet(m->variants[i], extra_snv, &temp, pair_cache);
        *out_result = res;
        if (res > 0)
        {
            can = false;
            break;
        }
    }

    mnv output = {};
    if (can)
    {
        output.variants.reserve(m->variants.size() + 1);
        for (int i = 0; i < m->variants.size(); i++)
        {
            output.variants.push_back(m->variants[i]);
        }
        output.variants.push_back(extra_snv);
        output.name = name_mnv(output.variants, true);

        if(blacklist_enable)
        {
            if(pair_blacklist.find(output.name) != pair_blacklist.end())
            {
                *out_result = 1;
                return output;
            }
        }

        *out_result = test_snv(output.variants, output.variants.size(), &output);
    }
    return output;
}

/*
    Parses a whole window of SNVs to doublets
    Every single combination of doublets is made
*/
mnv_window parse_window_to_doublets(snv_window& window, std::unordered_set<mnv, mnv_hash>& pair_cache, mnv_window& filtered, int window_id)
{
    std::vector<mnv> output_mnv;
    for(int i = 0; i < window.size(); i++)
    {
      for(int j = i + 1; j < window.size(); j++)
      {
        int result = 0;
        mnv m = make_doublet(window[i], window[j], &result, pair_cache);
        m.window_id = window_id;
        if(result == 0)
        {
          output_mnv.push_back(m);
        } else if(result == 2 && !settings.skip_filtered)
        {
            filtered.push_back(m);
        }
        pair_cache.insert(m);
      }
    }

    return output_mnv;
}

/*
    Parses a whole window of MNVs to the next phase size.
    If the specified size is larger than the amount of SNVs present in the window then this will not be performed.
*/
mnv_window parse_window_to_next_phase(mnv_window& mnv_w, snv_window& snv_w, std::unordered_set<mnv, mnv_hash>& pair_cache, mnv_window& filtered, int phase_size, int window_id)
{
  mnv_window output_mnv;

  std::unordered_set<std::string> name_cache;

  if(snv_w.size() < phase_size)
  {
    return output_mnv;
  }

  for(int i = 0; i < mnv_w.size(); i++)
  {
    for(int j = 0; j < snv_w.size(); j++)
    {
      if(!mnv_contains_snv(&mnv_w[i], snv_w[j]))
      {
        int result = 0;
        mnv m = make_next_phase(&mnv_w[i], snv_w[j], &result, pair_cache);
        if(result == 2 && !settings.skip_filtered)
        {
            filtered.push_back(m);
        }
        if(m.variants.size() > 0 && name_cache.find(m.name) == name_cache.end())
        {
            name_cache.insert(m.name);
            m.window_id = window_id;
            output_mnv.push_back(m);
        }
      }
    }
  }
  return output_mnv;
}

/*
    Makes all possible phase sizes from a window of SNVs, based on the specified size
    First we make all doublets, then starting from size = 3 we make every additional MNV until we reach either the SNV window size or the specified size.
*/
mnv_window make_all_phases_window(snv_window& window, std::unordered_set<mnv, mnv_hash>& pair_cache, mnv_window& filtered, int max_phase_size, int window_id)
{
  mnv_window dbl = parse_window_to_doublets(window, pair_cache, filtered, window_id);
  mnv_window cur_mnv = dbl;
  mnv_window result = dbl;
  for(int k = 3; k <= max_phase_size; k++)
  {
    mnv_window temp = parse_window_to_next_phase(cur_mnv, window, pair_cache, filtered, k, window_id);
    if(cur_mnv.empty())
    {
      break;
    }
    result.insert(result.end(), temp.begin(), temp.end());
    cur_mnv = temp;
  }
  
  return result;
}

/*
    Finds the actual exact position of a SNV in a read
    In the case of a deletion/insertion on the SNV position, the read is skipped
*/
int snv_relative_pos(bam1_t* b, int snv_pos)
{
    uint32_t* cigar = bam_get_cigar(b);
    int ref_pos = b->core.pos;
    int relative_pos = 0;

    for(uint32_t i = 0; i < b->core.n_cigar; ++i)
    {
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t oplen = bam_cigar_oplen(cigar[i]);
        uint32_t cons_ref = bam_cigar_type(op) & 2;
        uint32_t cons_query = bam_cigar_type(op) & 1;

        if(cons_ref && snv_pos >= ref_pos && snv_pos < ref_pos + oplen)
        {
            if(op == BAM_CDEL || op == BAM_CREF_SKIP)
            {
                return -1;
            }
            int offset = snv_pos - ref_pos;
            return relative_pos + offset;
        }

        if(cons_ref)
        {
            ref_pos += oplen;
        }

        if(cons_query)
        {
            relative_pos += oplen;
        }
    }
    return -1;
}

/*
    Loads reads for SNVs in a window
    Reads are stored as an unsigned integer hash based on read name
    Read names should be unique!
*/
void load_window_reads(samFile* bam, bam_hdr_t* bam_hdr, hts_idx_t* bam_idx, int chrom_id, snv_window& window)
{
    if(chrom_id < 0) return;
    
    std::hash<std::string> hasher;

    bam1_t* bam_read = bam_init1();
    int count = 0;
    int dup_counter = 0;
    for(snv* v : window)
    {
      int pos = v->pos;
      
      v->base_qual_sum = 0.0;
      v->loaded_reads = 1;
  
      hts_itr_t* iter = sam_itr_queryi(bam_idx, chrom_id, pos, pos + 1);
      if(!iter) return;
      while(sam_itr_next(bam, iter, bam_read) >= 0)
      {
        if(bam_read->core.qual < settings.min_read_quality) continue;
        
        int relative_pos = snv_relative_pos(bam_read, pos);
        if(relative_pos < 0 || relative_pos >= bam_read->core.l_qseq) continue;

        std::string read_name = bam_get_qname(bam_read);
        
        uint8_t* seq = bam_get_seq(bam_read);

        uint8_t* quals = bam_get_qual(bam_read);
        uint8_t q = quals[relative_pos];
        
        read r;
        r.read_name = read_name;
        r.quality = q;
        //unsigned int index = (unsigned int)hasher(read_name);

        v->covering_hashes.push_back(r);
        char base = seq_nt16_str[bam_seqi(seq, relative_pos)];

        if(base == v->alt)
        {
            if(q != 255)
            {
                v->base_qual_sum += ((float)(int)q / 10);
            }
            v->supporting_hashes.push_back(r);
        }
  
        count++;
      }

    //   if(v->base_qualities.size() > 1)
    //   {
    //     std::sort(v->base_qualities.begin(), v->base_qualities.end());
    //     int mid = v->base_qualities.size() / 2;
    //     mid = (mid % 2 == 0 ? mid : mid - 1);
    //     int qual = (int)v->base_qualities[mid];
    //     v->phred_qual = ((float)-qual / 10.0f); 

    //   } else
    //   {
    //     v->phred_qual = 0.0f;
    //   }

      if(v->supporting_hashes.size() > 0 && v->covering_hashes.size() > 0)
      {
        v->vaf = (float)v->supporting_hashes.size() / (float)v->covering_hashes.size();
      }

      v->mrd = v->supporting_hashes.size();
      v->dp = v->covering_hashes.size();
  
      bam_itr_destroy(iter);
    }
    bam_destroy1(bam_read);
}

/*
    Parses a whole SNV window into MNVs
*/
mnv_window parse_window(snv_window& window, mnv_window& filtered, std::string& in_bam, int chromosome, int window_id)
{
    mnv_window results;

    samFile* bam = sam_open(in_bam.c_str(), "r");
    if(bam == nullptr)
    {
      return results;
    }
    
    bam_hdr_t* bam_hdr = sam_hdr_read(bam);
    if(bam_hdr == nullptr)
    {
      return results;
    }
    
    hts_idx_t* bam_idx = sam_index_load(bam, in_bam.c_str());
    if(!bam_idx)
    {
      return results;
    }

    std::unordered_set<mnv, mnv_hash> cached_pairs;

    //Pre-allocate expected number of supporting/covering reads
    for(snv* v : window)
    {
        v->supporting_hashes.reserve(1024);
        v->covering_hashes.reserve(24576);
    }

    load_window_reads(bam, bam_hdr, bam_idx, chromosome, window);
    
    results = make_all_phases_window(window, cached_pairs, filtered, settings.mnv_size <= 1 ? window.size() : settings.mnv_size, window_id);

    for(snv* v : window)
    {
        if(v->loaded_reads && (!v->supporting_hashes.empty() || !v->covering_hashes.empty()))
        {
            v->supporting_hashes.clear();
            v->covering_hashes.clear();
            v->supporting_hashes = std::vector<read>(0);
            v->covering_hashes = std::vector<read>(0);
        }
    }
    
    sam_hdr_destroy(bam_hdr);
    hts_idx_destroy(bam_idx);
    sam_close(bam);

    std::stringstream strstr;
    strstr << "Finished window " << window_id << "(chr_id: " << chromosome << "). Total MNVs: " << results.size() << ", Total cached pairs: " << cached_pairs.size();
    log_info(strstr.str());

    return results;
}

/*
    Processes a single window on a single thread
    void* arg is used to pass in the chromosome and winow ID of this window
*/
void* process_window(void* arg)
{
    uint32_t packed_id = *(uint32_t*)arg;
    uint32_t chrom_id = (uint32_t) (packed_id >> 16);
    uint32_t window_id = (uint32_t) (packed_id & 0xFFFF); 

    if(chrom_id >= 0 && window_id >= 0 && window_id < total_windows)
    {
        all_mnv_windows[window_id] = parse_window(all_snv_windows[window_id], all_mnv_filtered[window_id], settings.in_bam, chrom_id, window_id);
    }

    return arg;
}

/*
    If specified, loads a blacklist (panel of normals) of MNVs that should be ignored.
*/
void load_blacklist(std::string file_path)
{
    log_info("Blacklist specified, trying to load it.");
    std::ifstream blacklist_file(file_path);
    if(!blacklist_file)
    {
        log_info("Blacklist file was specified, but failed to find/open the file. No blacklist applied.\n");
        return;
    }

    std::string line;
    while(std::getline(blacklist_file, line))
    {
        if(!line.empty())
        {

            if (line[line.length()-1] == '\n' || line[line.length()-1] == '\r') {
                line.erase(line.length()-1);
            }
            pair_blacklist.insert(line);
        }
    }

    blacklist_file.close();

    if(pair_blacklist.size() > 0)
    {
        for(auto& pair : pair_blacklist)
        {
            log_info("Found blacklist entry: " + pair);
        }
    }

    log_info("Found " + std::to_string(pair_blacklist.size()) + " pairs to be blacklisted.");
}

int main(int argc, char *argv[])
{
    settings = {};
    
    argparse::ArgumentParser program("MNV");

    program.add_argument("input_bam").help("Path to an input .bam file.").store_into(settings.in_bam);
    program.add_argument("input_vcf").help("Path to an input .vcf file.").store_into(settings.in_vcf);
    program.add_argument("output_dir").help("Path to the output file directory.").store_into(settings.out_path);
    program.add_argument("-O", "--output-name").default_value("results").help("The name for the output files. MNVs will be placed in <output_dir><name>.csv, the vcf in <output_dir><name>.vcf").store_into(settings.out_name);
    program.add_argument("-M", "--max-mnv-size").default_value(3).help("Determines the maximum amount of SNVs that can be in an MNV. Set to 0 for dynamic sizing based on window size, however this might increase runtime substantially.").store_into(settings.mnv_size);
    program.add_argument("-Q", "--read-quality").default_value(25).help("Only consider reads with an equal or higher mapping quality than the specified value.").store_into(settings.min_read_quality);
    program.add_argument("-V", "--verbose").default_value(false).help("Enables additional logging wile the program is running").store_into(settings.verbose);
    program.add_argument("-T", "--threads").default_value(4).help("Specifies the amount of threads the program should use. A higher number means an increase in work paralellization.").store_into(settings.num_threads);
    program.add_argument("-A", "--min-vaf-mnv").default_value(0.0001f).help("Minimum VAF for an MNV to be considered. MNVs with a lower VAF than this will not be output.").store_into(settings.min_vaf);
    program.add_argument("-S", "--min-vrd-snv").default_value(5).help("Minimum VRD for a SNV to be considered. SNVs with a lower VRD than this will be filtered prior to MNV analysis.").store_into(settings.min_mrd_snv);
    program.add_argument("-N", "--min-vrd-mnv").default_value(1).help("Minimum VRD for a MNV to be considered. MNVs with a lower VRD than this will not be output.").store_into(settings.min_mrd_mnv);
    program.add_argument("-Y", "--min-vaf-snv").default_value(0.0f).help("Minimum VAF for a SNV to be considered. SNVs with a lower VAF than this will be filtered prior to MNV analysis.").store_into(settings.min_snv_vaf);
    program.add_argument("-Z", "--max-vaf-snv").default_value(1.0f).help("Maximum VAF for a SNV to be considered. SNVs with a higher VAF than this will be filtered prior to MNV analysis.").store_into(settings.max_snv_vaf);
    program.add_argument("-B", "--min-bayesian").default_value(0.0).help("Minimum bayesian probability score for an MNV to be considered. MNVs with a lower bayesian probability than the specified value will be filtered.").store_into(settings.min_bayesian);
    program.add_argument("-F", "--min-phi").default_value(0.0).help("Minimum Phi-coefficient for a pair of SNVs to be considered. Pairs with a lower Phi-coefficient than the specified value will be discarded.").store_into(settings.min_phi);
    program.add_argument("-C", "--black-list").help("Path to a file containing a list of MNVs that should be ignored while making MNVs.").store_into(settings.blacklist_path);
    program.add_argument("-R", "--read-length").default_value(100).help("The maximum length in bp a read can be. SNVs will not be paired if their distance is larger than this value.").store_into(settings.read_length);
    program.add_argument("-J", "--min-jaccard").default_value(0.0).help("The minimum Jaccard index value for an MNV to be considered.").store_into(settings.jaccard);
    program.add_argument("-K", "--skip-filtered").default_value(false).help("Don't save filtered MNVs, only keep and output MNVs that passed the tests.").store_into(settings.skip_filtered);
    program.add_argument("-P", "--bayes-prior-mnv").default_value(0.5).help("Prior used for the Bayesian model. Change this to make SNV pairs less/more likely to be designated as real MNV. Default is 0.5 (no effect)").store_into(settings.bayes_prior);

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::exception &err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    try
    {
        std::filesystem::create_directories(settings.out_path);
        std::cout << "[MNVista] Output folder: " << settings.out_path << "\n";
        
        std::stringstream logfile;
        logfile << settings.out_path << std::filesystem::path::preferred_separator << settings.out_name << ".log";
        settings.log_stream.open(logfile.str());

        if(settings.log_stream.is_open());
        {
            std::cout << "[MNVista] Logging to:" << logfile.str() << std::endl;
        }
    }
    catch (std::exception &err)
    {
        std::cerr << err.what() << std::endl;
    }

    /*
        Print the used command in the logs
    */
    std::stringstream cmd_info;
    for(int i = 0; i < argc; i++)
    {
        cmd_info << argv[i] << " ";
    }
    log_info(cmd_info.str());


    if(settings.blacklist_path != "")
    {
        blacklist_enable = true;
        load_blacklist(settings.blacklist_path);
    }

    auto start = std::chrono::steady_clock::now();

    /*
        Start MNVista analysis
    */

    read_vcf(all_snvs, settings.in_vcf, settings.in_bam);

    total_windows = 0;

    for(int chr : chrom_contigs)
    {
        make_windows_chromosome(all_snvs, all_snv_windows, chr); 
    }

    total_windows = all_snv_windows.size();
    log_info("Created " + std::to_string(total_windows) + " total windows");
    
    all_mnv_windows = (mnv_window*)calloc(total_windows, sizeof(mnv_window));
    all_mnv_filtered = (mnv_window*)calloc(total_windows, sizeof(mnv_window));

    /*
        Set up multithreading based on input thread count setting
    */

    hts_tpool *p = hts_tpool_init(settings.num_threads);
    hts_tpool_process *q = hts_tpool_process_init(p, total_windows, 1);

    //Malloc array of packed chromosome/window IDs for each window
    uint32_t* packed_ids = (uint32_t*) calloc(total_windows, sizeof(uint32_t));

    log_info("Starting processing of windows...");

    for (int i = 0; i < total_windows; i++)
    {
        if (!all_snv_windows[i].empty())
        {
            //Pack the chromosome and window ID into a single unsigned integer to pass into the process_window function
            //Integer is unpacked within the process_window function again
            uint32_t chrom = (uint32_t)all_snv_windows[i][0]->chrom_id;
            uint32_t window_id = (uint32_t)i;
            packed_ids[i] = (chrom << 16) | (window_id & 0xFFFF);
            hts_tpool_dispatch(p, q, process_window, &packed_ids[i]);
        }
    }

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    //Write outputs
    write_mnv_list(all_mnv_windows, total_windows, false);
    write_mnv_list(all_mnv_filtered, total_windows, true);
    write_vcf_list(all_mnv_windows, vcf_contigs, total_windows);

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;

    std::chrono::duration<double> elapsed_seconds = end - start;
    double minutes = elapsed_seconds.count() / 60.0;
    log_info("Finished! Total execution time: " + std::to_string(minutes) + " minute(s).");

    free(all_mnv_windows);
    free(all_mnv_filtered);
    free(packed_ids);

    settings.log_stream.close();

    std::cout << "[MNVista] Finished run. \n";
    return 0;
}