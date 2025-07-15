#pragma once

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


#include "mnv.hpp"

float vaf_mean(const snv_window& variants, int num_snv);
float vaf_sd(const snv_window& variants, int num_snv, float* out_mean);

double test_phi(int num_both, int num_a, int num_b, int num_none);
double test_odds(int num_both, int num_a, int num_b, int num_none);
double test_bayesian(snv* snv_a, snv* snv_b, mnv* cur_mnv, int num_both, int num_a, int num_b, int num_none, double f_error, double f_haplo, double prior_mnv);