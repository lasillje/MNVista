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


#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm>

#include "utils.hpp"

#include "htslib/vcf.h"

/*
    Get the name of an MNV (chrN:pos1-pos2..:ref1>alt1-ref2>alt2...)
    If add_ref_alts is set to true then the reference/alt allele of each SNV is added to the name as well.
*/
std::string name_mnv(std::vector<snv *> &variants, bool add_ref_alts)
{
    if (!variants.empty())
    {
        std::stringstream strs;
        strs << variants[0]->chrom_name << ":";

        std::sort(variants.begin(), variants.end(), [](snv *a, snv *b)
                  { return a->pos < b->pos; });

        for (int i = 0; i < variants.size(); i++)
        {
            if (i > 0)
            {
                strs << "-";
            }
            strs << (variants[i]->pos + (settings.zero_indexed ? 0 : 1));
        }

        if (add_ref_alts)
        {
            strs << ":";
            for (int i = 0; i < variants.size(); i++)
            {
                if (i > 0)
                {
                    strs << "-";
                }
                strs << variants[i]->ref << ">" << variants[i]->alt;
            }
        }
        return strs.str();
    }
    return "Empty MNV";
}

/*
    Writes MNV list of filtered or unfiltered MNVs
*/
void write_mnv_list(mnv_window* windows, int total_windows, bool filtered)
{
    htsFile *vcf = bcf_open(settings.in_vcf.c_str(), "r");
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);

    std::stringstream file_name;
    file_name << settings.out_name;
    
    if(filtered)
    {
        file_name << "_filtered";
    }

    std::stringstream outstr;
    outstr << settings.out_path << std::filesystem::path::preferred_separator << file_name.str() << ".csv";
    std::string output_name = outstr.str();
    std::ofstream mnv_outfile;

    //Sort windows based on max MNV size
    //We want the biggest MNV on top for each window, and then decrease in size
    for(int i = 0; i < total_windows; i++)
    {
        if(!windows[i].empty())
        {
            std::sort(windows[i].begin(), windows[i].end(), [](const mnv &a, const mnv &b)
            {
                return a.variants.size() > b.variants.size(); });
        }
    }

    std::unordered_set<std::string> mnv_cache;

    mnv_outfile.open(output_name);
    if (mnv_outfile.is_open())
    {
        //Add header
        mnv_outfile << "WINDOW_ID;CHROM;MNV_NAME;VAF;SD;NUM_SUPPORTING;NUM_COVERING;PHI;COEFFICIENT_VARIATION;BAYESIAN;JACCARD;DISCORDANT;INDIVIDUAL_MUTATED;SIZE_MNV;DIST_MNV;QUALITIES;LOG_ODDS\n" << std::fixed;
        for(int i = 0; i < total_windows; i++)
        {
            if(windows[i].empty()) continue;
            for(int j = 0; j < windows[i].size(); j++)
            {
                for(const auto& m : windows[i])
                {   
                    if(mnv_cache.find(m.name) != mnv_cache.end() || m.variants.empty())
                    {
                        continue;
                    }
                    mnv_cache.insert(m.name);

                    int dist = m.variants[m.variants.size()-1]->pos - m.variants[0]->pos;
                    std::string chromosome = m.variants[0]->chrom_name;

                    std::stringstream discordants;
                    if(!m.discordant.empty())
                    {
                        for(int i = 0; i < m.variants.size(); i++)
                        {
                            discordants << m.discordant[i];
                            if(i != m.variants.size() - 1)
                            {
                                discordants << "|";
                            }
                        }
                    }

                    std::stringstream individual;
                    for(int i = 0; i < m.variants.size(); i++)
                    {
                        individual << m.variants[i]->mrd;
                        if(i != m.variants.size() - 1)
                        {
                            individual << "|";
                        }
                    }

                    std::stringstream qualities;
                    for(int i = 0; i < m.variants.size(); i++)
                    {
                        qualities << m.variants[i]->phred_qual * -10;
                        if(i != m.variants.size() - 1)
                        {
                            qualities << "|";
                        }
                    }
                    mnv_outfile << m.window_id << ";" << chromosome << ";" << m.name << ";" <<  m.vaf << ";" << m.sd << ";" << m.num_sup << ";" << m.num_cov << ";" << m.odds_phi << ";" << m.rsd << ";" << m.bayesian_prob << ";" << m.frac << ";" << discordants.str() << ";" << individual.str() << ";" << m.variants.size() << ";" << dist << ";" << qualities.str() << ";" << m.odds_ratio << "\n";
                }
            }
        }
    }
    mnv_outfile.close();
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
}

/*
    Write VCF based on available chrom contigs and unique SNVs
*/
void write_vcf_list(mnv_window* windows, const std::set<std::string>& contigs, int total_windows)
{
    htsFile *vcf = bcf_open(settings.in_vcf.c_str(), "r");
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);

    std::stringstream outstr;
    outstr << settings.out_path << std::filesystem::path::preferred_separator << settings.out_name << ".vcf";
    std::string output_name = outstr.str();
  
    std::ofstream vcf_outfile;

    std::set<snv*> unique_snv;
    for(int i = 0; i < total_windows; i++)
    {
        if(windows[i].empty());
        for(int j = 0; j < windows[i].size(); j++)
        {
            for(mnv& m : windows[i])
            {
                for(snv* v : m.variants)
                {
                    unique_snv.insert(v);
                }
            }
        }
    }

    std::vector<snv*> unique_variants_vec(unique_snv.begin(), unique_snv.end());

    std::sort(unique_variants_vec.begin(), unique_variants_vec.end(), [](snv* a, snv* b)
    {
      if(a->chrom_id != b->chrom_id)
      {
        return a->chrom_id < b->chrom_id;
      }
      return a->pos < b->pos;
    });
  
    vcf_outfile.open(output_name);
    if(vcf_outfile.is_open())
    {
        //Add header
        vcf_outfile << std::fixed;
        vcf_outfile << "##fileformat=VCFv4.2\n";
        vcf_outfile << "##source=MNV\n";
        vcf_outfile << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
        vcf_outfile << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n";
        vcf_outfile << "##INFO=<ID=MRD,Number=1,Type=Integer,Description=\"Mutant Read Depth\">\n";

        for(const auto& contig : contigs)
        {
            vcf_outfile << "##contig=<ID=" << contig << ">\n";
        }

        vcf_outfile << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
  
      for(snv* v : unique_variants_vec)
      {

        vcf_outfile << v->chrom_name << "\t" << v->pos + 1 << "\t" << "." << '\t' << v->ref << "\t" << v->alt << "\t" << "100" << "\tPASS" << "\t" << "AF=" << v->vaf << ";MRD=" << v->mrd << ";DP=" << v->dp << "\n";
      }
    } else
    {
      log_error("Unable to create file " + output_name + "\n");
    }
    vcf_outfile.close();

    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
}

/*
    Very simple helper functions for logging
*/

void log_info(std::string msg)
{
    settings.log_stream << "INFO: " << msg << std::endl;
}

void log_error(std::string msg)
{
    settings.log_stream << "ERROR: " << msg << std::endl;
}