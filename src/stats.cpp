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


#include <cmath>

#include <iostream>
#include <sstream>
#include "stats.hpp"
#include "utils.hpp"

float vaf_mean(const snv_window& variants, int num_snv)
{
    float vaf = 0;
    for (int i = 0; i < num_snv; i++)
    {
        vaf += variants[i]->vaf;
    }
    vaf /= (float)num_snv;
    return vaf;
}

/*
    Calculates the standard deviation and mean of the VAF from a list of variants.
    out_mean can be used as an additional output of the mean, if a float pointer is given.
*/
float vaf_sd(const snv_window& variants, int num_snv, float* out_mean)
{
    float mean = vaf_mean(variants, num_snv);

    if (out_mean != nullptr)
    {
        *out_mean = mean;
    }

    std::vector<float> devs;
    for (int i = 0; i < num_snv; i++)
    {
        float d = variants[i]->vaf - mean;
        devs.push_back(d * d);
    }

    float total = 0;
    for (int i = 0; i < num_snv; i++)
    {
        total += devs[i];
    }

    total /= (float)num_snv;
    return std::sqrt(total);
}

/*
    Log odds test, deprecated/unused
*/
double test_odds(int num_both, int num_a, int num_b, int num_none)
{
    double a = (double)(num_a + 1);
    double b = (double)(num_b + 1);
    double both = (double)(num_both + 1);
    double none = (double)(num_none + 1);

    double odds = (both * none) / (a * b);

    if(odds == 0)
    {
        return 0.0;
    }

    return std::log10(odds);
}

/*
    Phi coefficient test based on a 2x2 contingency table of read counts
*/
double test_phi(int num_both, int num_a, int num_b, int num_none)
{
    double total = num_both + num_a + num_b + num_none;

    long double alpha = (num_both + 1) / total;
    long double beta = (num_a + 1) / total;
    long double gamma = (num_b + 1) / total;
    long double delta = (num_none + 1) / total;

    long double numerator = (alpha * delta) - (beta * gamma);
    long double denominator = std::sqrt(((alpha + beta)*(alpha + gamma)*(beta + delta)*(gamma + delta)));

    if(settings.verbose)
    {
        std::stringstream debugstr;
        debugstr << "PHI A/B/C/D/NUM/DEM/RES: " << alpha << ", " << beta << ", " << gamma << ", " << delta << ", " << numerator << ", " << denominator << ", " << (numerator / denominator) << "\n";
        std::cerr << debugstr.str();
    }

    return static_cast<double>((numerator / denominator));
}

double test_bayesian(snv* snv_a, snv* snv_b, int num_both, int num_a, int num_b, int num_none, double f_error, double f_haplo, double prior_mnv, double prior_err)
{
    /*
        if the amount of reads containing both mutations is 0 we can stop early
        as the SNVs are mutually exclusive
    */
    if(num_both == 0)
    {
        return 0.0;
    }

    double A = (double)num_both; //Read counts for both alt allele
    double B = (double)num_a; //Read counts for only alt allele of snv_a
    double C = (double)num_b; //Read counts for only alt allele of snv_b
    double D = (double)num_none; //Read counts for both ref allele

    //Smoothing factor to prevent overflow/underflow
    //And to provide a smooth range of values.
    //The log likelihoods are calculated using raw read count integers to keep resolution
    //Maybe use long double instead of this?
    double a = 1.0 / (A + B + C);
    
    double eps_a = std::pow(10, snv_a->phred_qual); //Error rate for SNV a based on base quality. snv_a->phred_qual is equal to -Q / 10, in which Q is its base quality.
    double eps_b = std::pow(10, snv_b->phred_qual); //Same as above line but for snv_b

    double log_eps_a = std::log(eps_a);
    double log_eps_b = std::log(eps_b);

    //Error rate for A + B (all alt alleles for snv_a)
    double n_alt_a = A + B;
    double ll_a = n_alt_a * log_eps_a;

    //Error rate for A + C(all alt allels for snv_b)
    double n_alt_b = A + C;
    double ll_b = n_alt_b * log_eps_b;
    
    //Calculate chance for Model 1 (the two SNVs make a pair)
    double p_m1 =  A * std::log(f_haplo) + D * std::log(1.0 - f_haplo) + (B * log_eps_a + C * log_eps_b);

    //Calculate chance for Model 2A, Model 2B (Either one SNV is real and the other is caused by error)
    //Also consider D into this model in order to account for all available data
    double p_m2a = std::max((A + B) * std::log(f_error), ll_a) + ll_b + D * std::log(1-f_error);
    double p_m2b = std::max((A + C) * std::log(f_error), ll_b) + ll_a + D * std::log(1-f_error);
    double p_m2 = std::max(p_m2a, p_m2b);

    //Need to apply smoothing here otherwise the posterior collapses to either 0 or 1 due
    //to large read counts. Log likelihoods can reach up to -1000, which if put into the exponent
    //results in a number that cannot be put into a 64-bit floating point number.
    double lm1 = std::exp(p_m1 * a) * prior_mnv;
    double lm2 = std::exp(p_m2 * a) * (1.0 - prior_mnv);

    double denom = (lm1 + lm2);
    double post = lm1 / denom;

    if(settings.verbose)
    {
        std::stringstream debugstr;
        debugstr << snv_a->chrom_name << ":" << snv_a->pos << "-" << snv_b->pos << " -  A, B, C: " << A << ", " << B << ", " << C << " STATS: " << p_m1 << " , " << p_m2 << " (" << p_m2a << " , " << p_m2b << ")" << " ,  M1/M2/DENOM/POST: " << lm1 << " , " << lm2 << " , " << denom << " , " << post;
        log_info(debugstr.str());
    }

    return post;
}
