#pragma once

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


#include "mnv.hpp"

float vaf_mean(const snv_window& variants, int num_snv);
float vaf_sd(const snv_window& variants, int num_snv, float* out_mean);

double test_phi(int num_both, int num_a, int num_b, int num_none);
double test_odds(int num_both, int num_a, int num_b, int num_none);
double test_bayesian(snv* snv_a, snv* snv_b, int num_both, int num_a, int num_b, int num_none, double f_error, double f_haplo, double prior_mnv, double prior_err);