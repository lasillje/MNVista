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

std::string name_mnv(std::vector<snv*> &variants, bool add_ref_alts);

void write_mnv_list(mnv_window* windows, int total_windows, bool filtered);
void write_vcf_list(mnv_window* windows, const std::set<std::string>& contigs, int total_windows);

void log_info(std::string msg);
void log_error(std::string msg);