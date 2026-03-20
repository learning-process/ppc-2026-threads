#include "frolova_s_radix_sort_double/omp/include/ops_omp.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"

namespace frolova_s_radix_sort_double {

FrolovaSRadixSortDoubleOMP::FrolovaSRadixSortDoubleOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool FrolovaSRadixSortDoubleOMP::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleOMP::RunImpl() {
  const std::vector<double> &input = GetInput();
  if (input.empty()) {
    return false;
  }

  std::vector<double> working = input;

  const int radix = 256;      
  const int num_bits = 8;
  const int num_passes = sizeof(uint64_t);  

  for (int pass = 0; pass < num_passes; ++pass) {
    std::vector<int> count(radix, 0);

#pragma omp parallel
    {
      std::vector<int> local_count(radix, 0);
#pragma omp for nowait
      for (size_t i = 0; i < working.size(); ++i) {
        auto bits = std::bit_cast<uint64_t>(working[i]);
        int byte = static_cast<int>((bits >> (pass * num_bits)) & 0xFF);
        ++local_count[byte];
      }
#pragma omp critical
      {
        for (int j = 0; j < radix; ++j) {
          count[j] += local_count[j];
        }
      }
    }

    int total = 0;
    for (int i = 0; i < radix; ++i) {
      int old = count[i];
      count[i] = total;
      total += old;
    }

    std::vector<double> temp(working.size());
#pragma omp parallel for
    for (size_t i = 0; i < working.size(); ++i) {
      auto bits = std::bit_cast<uint64_t>(working[i]);
      int byte = static_cast<int>((bits >> (pass * num_bits)) & 0xFF);
      int pos = __sync_fetch_and_add(&count[byte], 1);
      temp[pos] = working[i];
    }

    working.swap(temp);
  }

  std::vector<double> negative;
  std::vector<double> positive;
  for (double val : working) {
    if (val < 0) {
      negative.push_back(val);
    } else {
      positive.push_back(val);
    }
  }
  std::ranges::reverse(negative);

  working.clear();
  working.insert(working.end(), negative.begin(), negative.end());
  working.insert(working.end(), positive.begin(), positive.end());

  GetOutput() = std::move(working);
  return true;
}

bool FrolovaSRadixSortDoubleOMP::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
