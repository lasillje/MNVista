#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm>

#include "utils.hpp"

#include "htslib/vcf.h"

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
    std::string outputName = outstr.str();
    std::ofstream mnv_outfile;

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

    mnv_outfile.open(outputName);
    if (mnv_outfile.is_open())
    {
        mnv_outfile << "WINDOW_ID;CHROM;MNV_NAME;VAF;SD;NUM_SUPPORTING;NUM_COVERING;PHI;COEFFICIENT_VARIATION;BAYESIAN;FRACTION_MUTATED;DISCORDANT;INDIVIDUAL_MUTATED;SIZE_MNV;DIST_MNV\n" << std::fixed;
        for(int i = 0; i < total_windows; i++)
        {
            if(windows[i].empty()) continue;
            for(int j = 0; j < windows[i].size(); j++)
            {
                bool has_outputted = false;
                for(const auto& m : windows[i])
                {   
                    if(mnv_cache.find(m.name) != mnv_cache.end() || m.variants.empty())
                    {
                        continue;
                    }
                    has_outputted = true;
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
                    mnv_outfile << m.window_id << ";" << chromosome << ";" << m.name << ";" <<  m.vaf << ";" << m.sd << ";" << m.num_sup << ";" << m.num_cov << ";" << m.odds_phi << ";" << m.rsd << ";" << m.bayesian_prob << ";" << m.frac << ";" << discordants.str() << ";" << individual.str() << ";" << m.variants.size() << ";" << dist << "\n";
                }
                if(has_outputted)
                {
                    mnv_outfile << "\n";
                }
            }
        }
    }
    mnv_outfile.close();
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
}

void write_vcf_list(mnv_window* windows, int total_windows)
{
    htsFile *vcf = bcf_open(settings.in_vcf.c_str(), "r");
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);

    std::stringstream outstr;
    outstr << settings.out_path << std::filesystem::path::preferred_separator << settings.out_name << ".vcf";
    std::string outputName = outstr.str();
  
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

    std::vector<snv*> uniqueVariantsVec(unique_snv.begin(), unique_snv.end());

    std::sort(uniqueVariantsVec.begin(), uniqueVariantsVec.end(), [](snv* a, snv* b)
    {
      if(a->chrom_id != b->chrom_id)
      {
        return a->chrom_id < b->chrom_id;
      }
      return a->pos < b->pos;
    });
  
    vcf_outfile.open(outputName);
    if(vcf_outfile.is_open())
    {
        vcf_outfile << std::fixed;
        vcf_outfile << "##fileformat=VCFv4.2\n";
        vcf_outfile << "##source=MNV\n";
        vcf_outfile << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
        vcf_outfile << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n";
        vcf_outfile << "##INFO=<ID=MRD,Number=1,Type=Integer,Description=\"Mutant Read Depth\">\n";
        vcf_outfile << "##contig=<ID=chr1>\n##contig=<ID=chr2>\n##contig=<ID=chr3>\n##contig=<ID=chr4>\n##contig=<ID=chr5>\n##contig=<ID=chr6>\n##contig=<ID=chr7>\n##contig=<ID=chr8>\n##contig=<ID=chr9>\n##contig=<ID=chr10>\n##contig=<ID=chr11>\n##contig=<ID=chr12>\n##contig=<ID=chr13>\n##contig=<ID=chr14>\n##contig=<ID=chr15>\n##contig=<ID=chr16>\n##contig=<ID=chr17>\n##contig=<ID=chr18>\n##contig=<ID=chr19>\n##contig=<ID=chr20>\n##contig=<ID=chr21>\n##contig=<ID=chr22>\n##contig=<ID=chrX>\n##contig=<ID=chrY>\n";
        vcf_outfile << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
  
      for(snv* v : uniqueVariantsVec)
      {

        vcf_outfile << v->chrom_name << "\t" << v->pos + 1 << "\t" << "." << '\t' << v->ref << "\t" << v->alt << "\t" << "100" << "\tPASS" << "\t" << "AF=" << v->vaf << ";MRD=" << v->mrd << ";DP=" << v->dp << "\n";
      }
  
  
    } else
    {
      log_error("Unable to create file " + outputName + "\n");
    }
    vcf_outfile.close();

    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
}

void log_info(std::string msg)
{
    settings.log_stream << "INFO: " << msg << std::endl;
}

void log_error(std::string msg)
{
    settings.log_stream << "ERROR: " << msg << std::endl;
}