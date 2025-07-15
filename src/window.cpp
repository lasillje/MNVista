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

#include "window.hpp"
#include "utils.hpp"

#include <iostream>
#include <algorithm>
#include <functional>

/*
    Makes windows of SNVs of within a continuous stretch of (default) 100bp of each other.
*/
void make_windows_chromosome(std::vector<snv> &variant_list, std::vector<snv_window>& output_windows, int chromosome)
{
    std::vector<snv_window> windows;
    snv_window cur_window;

    std::string cur_chrom_name;

    int last_pos = -1;
    
    for(int i = 0; i < variant_list.size(); i++)
    {
        int cur_chrom = variant_list[i].chrom_id;

        if(cur_chrom < chromosome) continue;
        if(cur_chrom > chromosome)
        {
            if(cur_window.size() > 1)
            {
                windows.push_back(cur_window);
                cur_window.clear();
            }
            break;
        }

        if(cur_chrom_name.empty())
        {
            cur_chrom_name = variant_list[i].chrom_name;
        }

        if(last_pos == -1)
        {
            cur_window.push_back(&variant_list[i]);
            last_pos = (int)variant_list[i].pos;
        } else
        {
            int dist = std::abs((int)last_pos - (int)variant_list[i].pos);
            if(dist < settings.read_length)
            {
                cur_window.push_back(&variant_list[i]);
                last_pos = (int)variant_list[i].pos;
            } else
            {
                if(cur_window.size() > 1)
                {
                    windows.push_back(cur_window);
                }
                cur_window.clear();
                cur_window.push_back(&variant_list[i]);
                last_pos = (int)variant_list[i].pos;
            }
        }
    }

    //For the very last window, if it is not empty then also add this
    if(!cur_window.empty())
    {
        windows.push_back(cur_window);
    }

    if(windows.empty())
    {
        return;
    } else
    {
        log_info("Found " + std::to_string(windows.size()) + " windows for " + cur_chrom_name);
    }

    std::sort(windows.begin(), windows.end(), [](const std::vector<snv *> &a, const std::vector<snv *> &b)
              {
                  int chromA = a[0]->chrom_id;
                  int chromB = b[0]->chrom_id;

                  if (chromA != chromB)
                  {
                      return chromA < chromB;
                  }
                  return a.size() > b.size(); });

    for(int i = 0; i < windows.size(); i++)
    {
        output_windows.push_back(windows[i]);
    }
}
