#pragma once

#include "mnv.hpp"

std::string name_mnv(std::vector<snv*> &variants, bool add_ref_alts);

void write_mnv_list(mnv_window* windows, int total_windows, bool filtered);
void write_vcf_list(mnv_window* windows, int total_windows);

void log_info(std::string msg);
void log_error(std::string msg);