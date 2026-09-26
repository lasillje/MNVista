#include "bamedit.hpp"

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/thread_pool.h"
#include "argparse/argparse.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <random>
#include <map>
#include <algorithm>
#include <mutex>

static std::mutex read_mtx;

std::string in_bed;
std::string in_ref;
std::string in_bam;

std::string out_path;
std::string out_name;

std::string out_vcf;
std::string out_bam;

static int MIN_MUTS = 3;
static int SEED = 1;
static int THREADS = 1;
static int CREATE_BAM = 1;
static bool VERBOSE = false;

static std::map<int, std::vector<snv*>> haplotypes;
static std::map<int, std::vector<std::string>> haplotype_reads;
static std::map<int, int> haplotype_chromosomes;

static std::map<std::string, std::vector<snv*>> reads_to_snvs;

std::vector<snv> load_spikes(const std::string& in_file)
{
    std::vector<snv> snvs;
    std::ifstream spike_file(in_file);

    if(!spike_file)
    {
        throw std::runtime_error("Could not open " + in_file);
    }

    std::string line;
    while(std::getline(spike_file, line))
    {
        if(line.empty() || line[0] == '#') continue;
        
        std::istringstream in_string(line);

        snv cur_snv;
        in_string >> cur_snv.chrom >> cur_snv.pos >> cur_snv.vaf >> cur_snv.haplotype;

        if(!in_string)
        {
            std::cerr << "WARNING: Skipped line " << line << " due to incorrect format.\n";
            continue;
        }
        cur_snv.vrd = -1;
        snvs.push_back(cur_snv);
    }

    std::sort(snvs.begin(), snvs.end(), [](snv& a, snv& b)
    {

      std::string subA = a.chrom.substr(3);
      std::string subB = b.chrom.substr(3);

      int chromA = subA == "X" ? 9998 : subA == "Y" ? 9999 : std::stoi(subA);
      int chromB = subB == "X" ? 9998 : subB == "Y" ? 9999 : std::stoi(subB);

      if(chromA != chromB)
      {
        return chromA < chromB;
      }
      return a.pos < b.pos;
    });


    return snvs;
}

void assign_bases(std::vector<snv>& snv_list, const std::string& ref_file)
{
    faidx_t* fai = fai_load(ref_file.c_str());
    int len = 0;

    std::mt19937 rng(SEED);
    std::uniform_int_distribution<int> dist(0, 2);

    for(snv& s : snv_list)
    {
        char* seq = faidx_fetch_seq(fai, s.chrom.c_str(), s.pos - 1, s.pos - 1, &len);
        
        if(len < 1)
        {
            std::cerr << "WARNING: Unable to find reference base for SNV " << s.chrom << " " << s.pos << ", skipping.\n";
            free(seq);
            continue;
        }

        char ref_base = seq[0];

        std::vector<char> choices = {'A','C','G','T'};
        std::remove(choices.begin(), choices.end(), ref_base);

        s.ref = ref_base;
        s.alt = choices[dist(rng)];

        free(seq);
    }

    fai_destroy(fai);
}

void assign_haplotypes(std::vector<snv>& snv_list)
{
    for(int i = 0; i < snv_list.size(); i++)
    {
        haplotypes[snv_list[i].haplotype].push_back(&snv_list[i]);
    }
}

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

void* assign_reads(void* arg)
{
    uint32_t hap = *(uint32_t*)arg;

    if(VERBOSE)
    {
        std::stringstream sstr;
        sstr << "Assigning reads for haplotype " << hap << std::endl;
        std::cout << sstr.str();
    }
    
    std::mt19937 rng(SEED);

    htsFile* in = sam_open(in_bam.c_str(), "r");

    if(in == nullptr)
    {   
        std::cerr << "Couldn't find input BAM file.";
        return arg;
    }

    hts_idx_t* idx = sam_index_load(in, in_bam.c_str());
    bam_hdr_t* hdr = sam_hdr_read(in);

    if(idx == nullptr || hdr == nullptr)
    {
        std::cerr << "Couldn't open index or header file.";
        return arg;
    }

    bam1_t* bam_read = bam_init1();

    int cur_haplotype = (int)hap;
    double cur_vaf = 0.0;

    if(haplotypes[hap].size() <= 0)
    {
        if(VERBOSE)
        {
            std::stringstream sstr;
            sstr << "Skipped assigning reads for haplotype " << hap << " due to no SNVs" << std::endl;
            std::cout << sstr.str();
        }
        return arg;
    }

    std::vector<std::set<std::string>> read_sets(haplotypes[hap].size());
    int counter = 0;

    for(snv* s : haplotypes[hap])
    {
        if(s == nullptr)
        {
            continue;
        }

        s->tid = bam_name2id(hdr, s->chrom.c_str());
        cur_vaf = s->vaf;

        int dp = 0;

        hts_itr_t* iter = sam_itr_queryi(idx, s->tid, s->pos - 1, s->pos);
        if(!iter) continue;
        while(sam_itr_next(in, iter, bam_read) >= 0)
        {

            int relative_pos = snv_relative_pos(bam_read, s->pos - 1);
            if(relative_pos < 0 || relative_pos >= bam_read->core.l_qseq) continue;
        
            std::string read_name = bam_get_qname(bam_read);

            read_sets[counter].insert(read_name);
            dp++;
        }

        s->dp = dp;

        if(counter == 0)
        {
            haplotype_chromosomes[s->haplotype] = s->tid;
        }

        bam_itr_destroy(iter);
        counter++;
    }

    std::set<std::string> intersect_set(read_sets[0]);
    
    for (int i = 1; i < read_sets.size(); i++)
    {
        std::set<std::string> temp;
        std::set_intersection(intersect_set.begin(), intersect_set.end(), read_sets[i].begin(), read_sets[i].end(), std::inserter(temp, temp.begin()));
        intersect_set = temp;
    }

    std::vector<std::string> unique_reads(intersect_set.begin(), intersect_set.end());

    if(unique_reads.size() <= 0)
    {
        std::stringstream sstr;
        sstr << "WARNING: No common reads found for haplotype " << hap << ", no spikes can be made." << std::endl;
        std::cout << sstr.str();
    }

    std::vector<std::vector<std::string>> mutated_reads_per_snv;

    std::vector<snv*> temp_snvs;
    for(int i = 0; i < haplotypes[hap].size(); i++)
    {
        temp_snvs.push_back(haplotypes[hap][i]);
    }

    
    std::sort(temp_snvs.begin(), temp_snvs.end(), [](snv* a, snv* b)
    {
        return a->vaf > b->vaf;
    });

    
    for(int i = 0; i < temp_snvs.size(); i++)
    {
        int num_mutated = MIN_MUTS;

        if(i == 0)
        {
            num_mutated = std::round(temp_snvs[i]->vaf * unique_reads.size());
        } else
        {
            num_mutated = std::round((temp_snvs[i]->vaf / temp_snvs[i-1]->vaf) * mutated_reads_per_snv[i-1].size());
        } 

        if(num_mutated <  MIN_MUTS)
        {
            num_mutated = MIN_MUTS;
        }

        std::vector<std::string> mutated_reads;
        std::ranges::sample(i == 0 ? unique_reads : mutated_reads_per_snv[i-1], std::back_inserter(mutated_reads), num_mutated, rng);
        mutated_reads_per_snv.push_back(mutated_reads);
    }

    //haplotype_reads[cur_haplotype] = mutated_reads;

    {
        std::lock_guard<std::mutex> lock(read_mtx);

        int count = 0;
        for(snv* s : haplotypes[hap])
        {
            for(int i = 0; i < mutated_reads_per_snv[count].size(); i++)
            {
                std::string readname = mutated_reads_per_snv[count][i];

                if(reads_to_snvs.find(readname) == reads_to_snvs.end())
                {
                    reads_to_snvs[readname] = {};
                }

                reads_to_snvs[readname].push_back(s);
            }
            count++;
        }
    }

    int count = 0;
    for(snv* s : haplotypes[hap])
    {
        if(s == nullptr)
        {
            continue;
        }
        s->vrd = 0;//mutated_reads_per_snv[count].size();
        count++;
    }
    
    if(VERBOSE)
    {
        std::stringstream sstr;
        sstr << "Finished haplotype " << hap << std::endl;
        std::cout << sstr.str();
    }

    bam_destroy1(bam_read);
    sam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(in);

    return arg;
}

void spike_and_write_reads(const std::string& in_bam, const std::string& out_bam)
{
    htsFile* in = sam_open(in_bam.c_str(), "r");
    if(in == nullptr)
    {
        std::cerr << "Couldn't find input BAM file." << std::endl;
        return;
    }
    hts_idx_t* idx = sam_index_load(in, in_bam.c_str());
    bam_hdr_t* hdr = sam_hdr_read(in);

    htsFile *out = sam_open(out_bam.c_str(), "wb");
    int hdr_res = sam_hdr_write(out, hdr);

    if(out == nullptr)
    {
        std::cerr << "Couldn't open output BAM file." << std::endl;
        return;
    }

    if(hdr_res < 0)
    {
        std::cerr << "WARNING: Failed to write SAM header." << std::endl;
    }

    bam1_t* b = bam_init1();
    while(sam_read1(in, hdr, b) >= 0)
    {
        std::string qname = bam_get_qname(b);

        if(reads_to_snvs.find(qname) != reads_to_snvs.end())
        {
            for(snv* s : reads_to_snvs[qname])
            {
                if(s->tid == b->core.tid)
                {
                    int qpos = snv_relative_pos(b, s->pos - 1);
                    if(qpos < 0 || qpos >= b->core.l_qseq) continue;

                    uint8_t* seq = bam_get_seq(b);
                    uint8_t nt16 = seq_nt16_table[(int)s->alt];
                    bam_set_seqi(seq, qpos, nt16);
                    s->vrd += 1;
                }
            }
        }

        // std::vector<int> hap_containing;

        // for(int i = 0; i < haplotype_reads.size(); i++)
        // {   
        //     if(haplotype_chromosomes[i] ==  b->core.tid)
        //     {
        //         if(std::find(haplotype_reads[i].begin(), haplotype_reads[i].end(), qname) != haplotype_reads[i].end())
        //         {
        //             hap_containing.push_back(i);
        //         }
        //     }
        // }

        // if(hap_containing.size() >= 1)
        // {
        //     for(int i = 0; i < hap_containing.size(); i++)
        //     {
        //         for(snv* s : haplotypes[hap_containing[i]])
        //         {
        //             if(s->tid == b->core.tid)
        //             {
        //                 int qpos = snv_relative_pos(b, s->pos - 1);
        //                 if(qpos < 0 || qpos >= b->core.l_qseq) continue;

        //                 uint8_t* seq = bam_get_seq(b);
        //                 uint8_t nt16 = seq_nt16_table[(int)s->alt];
        //                 bam_set_seqi(seq, qpos, nt16);
        //             }
        //         }
        //     }
        // }
        int res = sam_write1(out, hdr, b);
        if(res < 0)
        {
            std::cerr << "WARNING: Failed to write SAM record for read " << qname << "\n";
        }
    }

    sam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(in);

    sam_close(out);
}

void write_vcf(std::vector<snv>& snv_list)
{
    std::ofstream vcf_file(out_vcf);

    if(!vcf_file)
    {
        throw std::runtime_error("Could not open " + out_vcf);
    }

    vcf_file << std::fixed;
    vcf_file << "##fileformat=VCFv4.2\n";
    vcf_file << "##source=BamEdit\n";
    vcf_file << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
    vcf_file << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n";
    vcf_file << "##INFO=<ID=VRD,Number=1,Type=Integer,Description=\"Variant Read Depth\">\n";
    vcf_file << "##INFO=<ID=HAPLO,Number=1,Type=Integer,Description=\"Haplotype identifier used by BamEdit\">\n";

    std::set<std::string> unique_chroms;
    for(int i = 0; i < haplotypes.size(); i++)
    {
        if(haplotypes[i].size() > 0)
        {   
            unique_chroms.insert(haplotypes[i][0]->chrom);
        }
    }


    for(const auto& chrom : unique_chroms)
    {
        vcf_file << "##contig=<ID=" << chrom << ">\n";
    }

    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    for(snv& s : snv_list)
    {
        if(s.vrd < 0)
        {
            std::cerr << "Skipped VCF output of " << s.chrom << " " << s.pos << " due to VRD of " << s.vrd << std::endl;
            continue;
        }
        vcf_file << s.chrom << "\t" << s.pos << "\t" << "." << "\t" << s.ref << "\t" << s.alt << "\t" << "." << "\t" << "PASS" << "\t" << "AF=" << s.vaf << ";VRD=" << s.vrd << ";DP=" << s.dp << ";HAPLO=" << s.haplotype << "\n";
    }

}

int main(int argc, char* argv[])
{
    argparse::ArgumentParser program("BamEdit");

    program.add_argument("input_bed").help("Path to an input .bed file.").store_into(in_bed);
    program.add_argument("input_ref").help("Path to an input indexed fasta file.").store_into(in_ref);
    program.add_argument("input_bam").help("Path to an input indexed bam file.").store_into(in_bam);
    program.add_argument("output_dir").help("Path to the output file directory.").store_into(out_path);
    program.add_argument("-O", "--output-name").default_value("results").help("The name for the output files. A newly constructed BAM will be placed in <output_dir>/<output_name>.bam. The output VCF will be placed in the same directory with the same name.").store_into(out_name);
    program.add_argument("-M", "--minimum-mutations").default_value(3).help("The minimum number of reads to mutate, incase the VAF is too low.").store_into(MIN_MUTS);
    program.add_argument("-S", "--seed").default_value(1).help("Sets the seed for random number generation. Keep this the same across runs for reproducability.").store_into(SEED);
    program.add_argument("-T", "--threads").default_value(1).help("The amount of threads the program uses for assigning reads to each SNV.").store_into(THREADS);
    program.add_argument("-B", "--output-bam").default_value(1).help("Set to 0 if you want to skip the creation of a new BAM file containing the mutated reads. Only the VCF with mutated reads will be output.").store_into(CREATE_BAM);
    program.add_argument("-V", "--verbose").default_value(false).help("Enable for more verbose logging.").store_into(VERBOSE);

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

    std::stringstream vcf_stream;
    vcf_stream << out_path << std::filesystem::path::preferred_separator << out_name << ".vcf";
    out_vcf = vcf_stream.str();

    std::stringstream bam_stream;
    bam_stream << out_path << std::filesystem::path::preferred_separator << out_name << ".bam";
    out_bam = bam_stream.str();

    std::cout << "Loading SNVs...\n";
    std::vector<snv> spikes = load_spikes(in_bed);
    std::cout << "Found " << spikes.size() << " SNVs" << std::endl;

    std::cout << "Assigning bases to SNVs..." << std::endl;
    assign_bases(spikes, in_ref);

    std::cout << "Assigning haplotypes to SNVs..." << std::endl;
    assign_haplotypes(spikes);
    std::cout << "Found " << haplotypes.size() << " haplotypes" << std::endl;

    std::cout << "Assigning reads to SNVs: " << THREADS << " threads to use." << std::endl;
    hts_tpool *p = hts_tpool_init(THREADS);
    hts_tpool_process *q = hts_tpool_process_init(p, haplotypes.size(), 1);

    uint32_t* haplotype_ids = (uint32_t*)calloc(haplotypes.size(), sizeof(uint32_t));

    for(int i = 0; i < haplotypes.size(); i++)
    {
        haplotype_ids[i] = i;
        hts_tpool_dispatch(p, q, assign_reads, &haplotype_ids[i]);
    }

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    if(CREATE_BAM == 1)
    {
        std::cout << "Writing new BAM file with modified reads..." << std::endl;
        spike_and_write_reads(in_bam, out_bam);
    }

    std::cout << "Writing VCF..." << std::endl;
    write_vcf(spikes);

    free(haplotype_ids);

    return 0;
}