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

#include <cmath>

#include <iostream>
#include <sstream>
#include "stats.hpp"
#include "utils.hpp"

float vaf_mean(const snv_window& variants, int num_snv)
{
    float vaf = 0;
    for(int i = 0; i < num_snv; i++)
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

    if(out_mean != nullptr)
    {
        *out_mean = mean;
    }

    std::vector<float> devs;
    for(int i = 0; i < num_snv; i++)
    {
        float d = variants[i]->vaf - mean;
        devs.push_back(d * d);
    }

    float total = 0;
    for(int i = 0; i < num_snv; i++)
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
    long double denominator = std::sqrt(((alpha + beta) * (alpha + gamma) * (beta + delta) * (gamma + delta)));

    return static_cast<double>((numerator / denominator));
}

//log10 gamma (Stirling approx)
double lgamma10(double z)
{
    return std::lgamma(z) / M_LN10;
}

//Log sum exp helper for numerical stability
double log10_sum_exp(double A, double B)
{
    double M = std::max(A, B);

    double eA = std::pow(10.0, A - M);
    double eB = std::pow(10.0, B - M);

    return M + std::log10(eA + eB);
}


double test_bayesian(mnv* cur_mnv, int num_both, int num_alt_1, int num_alt_2, int num_none,
                     double p_err, double prior_mnv)
{
	// skip calculation if there are no reads in A anyways, or if somehow MNVs of size >= 2 were entered
    if(num_both == 0 || cur_mnv->qualities.size() < 2 || cur_mnv->discordant_qualities.size() < 2)
    {
        return 0.0;
    }

    const double A = (double)num_both; // alt at both positions
    const double B = (double)num_alt_1; // alt at position 1 only
    const double C = (double)num_alt_2; // alt at position 2 only
    const double D = (double)num_none; // ref at both positions
    const double R = A + B + C + D;

    // accumulated Phred error scores
    const double E_A1 = cur_mnv->qualities[0];
    const double E_A2 = cur_mnv->qualities[1];
    const double E_B1 = cur_mnv->discordant_qualities[0];
    const double E_C2 = cur_mnv->discordant_qualities[1];

    //Dirichlet(1,1)
    const double fa = 1.0, fb = 1.0;
    const double logB0 = lgamma10(fa) + lgamma10(fb) - lgamma10(fa + fb);

    // Model 1: both alts lie on same MNV haplotype
    // L_M1 = 10^(-E_B1 - E_C2) * Beta(A + fa, R - A + fb) / Beta(fa, fb)
	// Here in log space for numerical stability
    const double log_L_M1 = -E_B1 - E_C2
                          + lgamma10(A + fa) + lgamma10(R - A + fb) - lgamma10(R + fa + fb)
                          - logB0;

    // Model 2: both positions are independent
	// numerical sanity checks for p_err
    const double log_p_err = (p_err <= 0.0)  ? -1.0e300 : std::log10(p_err);
    const double log_p_snv = (p_err >= 1.0)  ? -1.0e300 : std::log10(1.0 - p_err);

    const double log_L_ERR1 = -E_A1 - E_B1;
    const double log_L_ERR2 = -E_A2 - E_C2;

    const double log_L_SNV1 = lgamma10(A + B + fa) + lgamma10(C + D + fb)
                            - lgamma10(R + fa + fb) - logB0;
    const double log_L_SNV2 = lgamma10(A + C + fa) + lgamma10(B + D + fb)
                            - lgamma10(R + fa + fb) - logB0;

    const double log_L_M2 = log10_sum_exp(log_p_err + log_L_ERR1, log_p_snv + log_L_SNV1)
                          + log10_sum_exp(log_p_err + log_L_ERR2, log_p_snv + log_L_SNV2);

    // Posterior also in log space
    const double diff = (log_L_M2 + std::log10(1.0 - prior_mnv))
                      - (log_L_M1 + std::log10(prior_mnv));
	
	// Max value clamp for numerical stability/ prevent overflow 
    if(diff >  300.0) return 0.0;
    if(diff < -300.0) return 1.0;
    return 1.0 / (1.0 + std::pow(10.0, diff));
}
