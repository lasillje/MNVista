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

    return static_cast<double>((numerator / denominator));
}

double softmax(double A, double B)
{
    double M = std::max(A, B);

    double eA = std::pow(10.0, A - M);
    double eB = std::pow(10.0, B - M);

    return M + std::log10(eA + eB);
}


double test_bayesian(snv* snv_a, snv* snv_b, mnv* cur_mnv, int num_both, int num_a, int num_b, int num_none, double f_error, double f_haplo, double prior_mnv)
{
    if(num_both == 0 || cur_mnv->qualities.size() < 2)
    {
        return 0.0;
    }

    double A = (double)num_both; //Read counts for both alt allele
    double B = (double)num_a; //Read counts for only alt allele of snv_a
    double C = (double)num_b; //Read counts for only alt allele of snv_b
    double D = (double)num_none; //Read counts for both ref allele

    double e_b = cur_mnv->discordant_qualities[0];
    double e_c = cur_mnv->discordant_qualities[1];
    
    double e_a1 = cur_mnv->qualities[0];
    double e_a2 = cur_mnv->qualities[1];

    double T = A + B + C + D;
    
    double fa = A / T;

    double fa1 = (A + B + 1) / T;
    double fa2 = (A + C + 1) / T;
    
    double fr1 = (C + D + 1) / T;
    double fr2 = (B + D + 1) / T;

    double p_m1 = A * std::log10(fa) + D * std::log10(1.0 - fa) - e_b - e_c;
    
    double p_m2b = softmax(-e_a1 - e_b, ((A+B) * std::log10(fa1)) + ((C+D) * std::log10(fr1)));
    double p_m2c = softmax(-e_a2 - e_c, ((A+C) * std::log10(fa2)) + ((B+D) * std::log10(fr2)));

    double p_m2 = p_m2b + p_m2c;

    double log_prior = std::log10(prior_mnv);
    double log_not_prior = std::log10(1.0 - prior_mnv);

    double f_p_m1 = p_m1 + log_prior;
    double f_p_m2 = p_m2 + log_not_prior;
    double diff = f_p_m2 - f_p_m1;
    
    double post = 0.0;

    if(diff > 100.0)
    {
        post = 0.0;
    } else if (diff < -100.0)
    {
        post = 1.0;
    } else
    {
        post = 1.0 / (1.0 + std::pow(10.0, diff));
    }

    return post;
}
