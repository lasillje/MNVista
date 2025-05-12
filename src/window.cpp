#include "window.hpp"
#include "utils.hpp"

#include <iostream>
#include <algorithm>
#include <functional>

static bool compare_snv_windows(snv_window& a, snv_window& b)
{
    std::set<int> setA;
    std::set<int> setB;
  
    for(int i = 0; i < a.size(); i++)
    {
      setA.insert(a[i]->pos);
    }
  
    for(int i = 0; i < b.size(); i++)
    {
      setB.insert(b[i]->pos);
    }
  
    if(std::includes(setB.begin(), setB.end(), setA.begin(), setA.end()))
    {
      return true;
    }
    return false;
}

static float similarity_snv_windows(snv_window& a, snv_window& b)
{
    auto a_sorted = a;
    auto b_sorted = b;
  
    std::vector<snv*> intersect;
    std::set_intersection(a_sorted.begin(), a_sorted.end(), b_sorted.begin(), b_sorted.end(), std::back_inserter(intersect));
  
    std::vector<snv*> uni;
    std::set_union(a_sorted.begin(), a_sorted.end(), b_sorted.begin(), b_sorted.end(), std::back_inserter(uni));
  
    if(uni.empty())
    {
      return 1.0f;
    }
    float res = (float) intersect.size() / (float)uni.size();

    return res;
}

static std::vector<snv_window> filter_windows(std::vector<snv_window>& windows)
{
    std::vector<snv_window> filtered_windows;

    for (int i = 0; i < windows.size(); i++)
    {
        bool canAdd = true;
        for (int j = i; j >= 0; j--)
        {
            if (i == j)
                continue;
            if (windows[i][0]->chrom_id == windows[j][0]->chrom_id)
            {
                if (compare_snv_windows(windows[i], windows[j]))
                {
                    canAdd = false;
                    break;
                }
            }
        }
        if (canAdd)
        {
            filtered_windows.push_back(windows[i]);
        }
    }
    return filtered_windows;
}

static std::vector<snv_window> merge_windows(std::vector<snv_window>& windows)
{
    std::vector<snv_window> merged_windows;
    std::set<int> mergedIndices;
    for (int i = 0; i < windows.size(); i++)
    {
        if (mergedIndices.find(i) != mergedIndices.end())
            continue;
        std::sort(windows[i].begin(), windows[i].end());
        std::set<snv *> localUniqueSNV(windows[i].begin(), windows[i].end());
        for (int j = i + 1; j < windows.size(); j++)
        {
            if (mergedIndices.find(j) != mergedIndices.end())
                continue;
            if (windows[i][0]->chrom_id != windows[j][0]->chrom_id)
                break;
            std::sort(windows[j].begin(), windows[j].end());
            if (similarity_snv_windows(windows[i], windows[j]) > 0.25f)
            {
                // i = j;
                for(int x = 0; x < windows[j].size(); x++)
                {
                    localUniqueSNV.insert(windows[j][x]);
                }
                //localUniqueSNV.insert(filteredWindows[j].begin(), filteredWindows[j].end());
                mergedIndices.insert(j);
            }
        }
        snv_window finalVec(localUniqueSNV.begin(), localUniqueSNV.end());
        merged_windows.push_back(finalVec);
    }
    return merged_windows;
}

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

    log_info("Found " + std::to_string(windows.size()) + " windows for " + cur_chrom_name + "\n");

    if(windows.empty())
    {
        return;
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

void make_windows_chromosome_old(std::vector<snv> &variant_list, std::vector<snv_window>& output_windows, int chromosome)
{
    std::vector<snv_window> windows;
    std::unordered_set<std::string> unique_window_names;

    for (int i = 0; i < variant_list.size(); i++)
    {
        int curChromID = variant_list[i].chrom_id;

        if (curChromID > chromosome)
            break;
        if (curChromID < chromosome)
            continue;

        std::vector<snv *> curWindow;
        curWindow.push_back(&variant_list[i]);

        for (int j = 0; j < variant_list.size(); j++)
        {
            if (variant_list[j].chrom_id < curChromID || i == j)
                continue;

            bool isSame = true;
            int dist = std::abs((int)variant_list[i].pos - (int)variant_list[j].pos);
            if(dist == 0)
            {
                continue;
            }
            if (curChromID == variant_list[j].chrom_id)
            {
                if (dist < settings.read_length)
                {
                    curWindow.push_back(&variant_list[j]);

                }
                else if (j > i)
                {
                    isSame = false;
                }
            }
            else
            {
                isSame = false;
            }

            if (!isSame || j == variant_list.size() - 1)
            {
                if (curWindow.size() > 1)
                {
                    std::string windowName = name_mnv(curWindow, true);
                    if (unique_window_names.find(windowName) == unique_window_names.end())
                    {
                        unique_window_names.insert(windowName);
                        windows.push_back(curWindow);
                    }
                }
                break;
            }
        }
    }

    if(windows.empty())
    {
        return;
    }

    unique_window_names.clear();

    std::sort(windows.begin(), windows.end(), [](const std::vector<snv *> &a, const std::vector<snv *> &b)
              {
                  int chromA = a[0]->chrom_id;
                  int chromB = b[0]->chrom_id;

                  if (chromA != chromB)
                  {
                      return chromA < chromB;
                  }
                  return a.size() > b.size(); });

    std::vector<snv_window> filtered_windows = filter_windows(windows);

    std::cout << "Kept " << filtered_windows.size() << " windows after pruning.\n";

    windows.clear();
    
    std::vector<snv_window> merged_windows = merge_windows(filtered_windows);
    filtered_windows.clear();
    
    std::cout << "Kept " << merged_windows.size() << " windows after merging.\n";

    std::sort(merged_windows.begin(), merged_windows.end(), [](const std::vector<snv *> &a, const std::vector<snv *> &b)
    {
        int chromA = a[0]->chrom_id;
        int chromB = b[0]->chrom_id;
  
        if(chromA != chromB)
        {
            return chromA < chromB;
        }

        return a.size() > b.size();
     });

      for(int i = 0; i < merged_windows.size(); i++)
      {
        output_windows.push_back(merged_windows[i]);
      }

}