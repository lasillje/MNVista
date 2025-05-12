#pragma once

#include "mnv.hpp"

float vaf_mean(const snv_window& variants, int num_snv);
float vaf_sd(const snv_window& variants, int num_snv, float* out_mean);

double test_phi(int num_both, int num_a, int num_b, int num_none);
double test_odds(int num_both, int num_a, int num_b, int num_none);
double test_bayesian(snv* snv_a, snv* snv_b, int num_both, int num_a, int num_b, int num_none, double rsd, double frac_mutated);