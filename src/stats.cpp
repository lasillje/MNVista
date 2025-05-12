#include <cmath>

#include "stats.hpp"


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

    double odds = (double)(((double)num_both) * ((double)num_none)) / (double)(a * b);

    return std::log10(odds);
}

/*
    Phi coefficient test based on a 2x2 contingency table of read counts
*/
double test_phi(int num_both, int num_a, int num_b, int num_none)
{
    double total = num_both + num_a + num_b + num_none;

    double alpha = (num_both + 1) / total;
    double beta = (num_a + 1) / total;
    double gamma = (num_b + 1) / total;
    double delta = (num_none) / total;

    double numerator = (alpha * delta) - (beta * gamma);
    double denominator = std::sqrt((double)((alpha + beta)*(alpha + gamma)*(beta + delta)*(gamma + delta)));
    return numerator / denominator;
}

/*
    Calculates Bayesian probability for a given SNV pair based on various statistics
*/
double test_bayesian(snv* snv_a, snv* snv_b, int num_both, int num_a, int num_b, int num_none, double rsd, double frac_mutated)
{
    /*
        M1: H0 = ref, ref   [1 - f]
            H1 = alt, alt   [f]

        M2: H0 = ref, ref   [1 - f1 - f2]
            H1 = ref, alt   [f1]
            H2 = alt, ref   [f2]
        

        M1: data is in-phase
        M2: data is not in phase
        f: frequency of hypothesis, assumed to be the VAF of each mutation
           for the paired SNV the vaf is assumed to be the number of reads containing both mutations divided by the total number of overlapping reads.

        Here we assign probabilities to M1 and M2 as the sum of the probabilities of their hypotheses.
        For M1, it is assumed that the only explanation for SNVs is that they are in phase.
        For M2, it is assumed that the only explanation for SNVs is that they are present due to random error.

        Furthermore, we have to subtract the error rate from the probability of M1. This is calculated as:
        (B * -Q/10) + (C * -Q/10)

        In which Q is the median phred quality score per SNV.

    */

    double prior_in_phase = (1.0 - rsd); //* frac_mutated;

    double A = (double)num_both;
    double B = (double)num_a;
    double C = (double)num_b;
    double D = (double)num_none;
    double T = A + B + C + D + 2.0;

    double m1_f_a = ((A + 1) / T);
    double m1_f_d = 1.0 - m1_f_a;
    
    double m2_f_b = ((B + 1) / T);
    double m2_f_c = ((C + 1) / T);
    double m2_f_d = 1.0 - m2_f_b - m2_f_c;

    double p_m1 = (D * std::log10(m1_f_d)) + (A * std::log10(m1_f_a)) - (B * snv_a->phred_qual + C * snv_b->phred_qual);
    double p_m2 = (D * std::log10(m2_f_d)) + (B * std::log10(m2_f_b)) + (C * std::log10(m2_f_c));

    double numerator = p_m1 * prior_in_phase;

    double denominator = numerator + p_m2 * (1.0 - prior_in_phase);

    return numerator / denominator;
}