#include "chernov_t_radix_sort/tbb/include/ops_tbb.hpp"

#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"

namespace chernov_t_radix_sort {

ChernovTRadixSortTBB::ChernovTRadixSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ChernovTRadixSortTBB::ValidationImpl() {
  return true;
}

bool ChernovTRadixSortTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

constexpr int kBitsPerDigit = 8;
constexpr int kRadix = 1 << kBitsPerDigit;
constexpr uint32_t kSignMask = 0x80000000U;

void ChernovTRadixSortTBB::ComputePrefixSums(std::vector<int> &count) {
  for (size_t idx = 1; idx < count.size(); ++idx) {
    count[idx] += count[idx - 1];
  }
}

void ChernovTRadixSortTBB::RadixSortLSD(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  const size_t n = data.size();
  std::vector<uint32_t> temp(n);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&temp, &data](const tbb::blocked_range<size_t> &r) {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      temp[i] = static_cast<uint32_t>(data[i]) ^ kSignMask;
    }
  });

  std::vector<uint32_t> buffer(n);

  for (int byte_index = 0; byte_index < 4; ++byte_index) {
    const int shift = byte_index * kBitsPerDigit;

    tbb::combinable<std::vector<int>> local_counts([]() { return std::vector<int>(kRadix, 0); });

    tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
                      [&local_counts, &temp, shift](const tbb::blocked_range<size_t> &r) {
      auto &my_count = local_counts.local();
      for (size_t i = r.begin(); i != r.end(); ++i) {
        int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
        ++my_count[static_cast<size_t>(digit)];
      }
    });

    std::vector<int> count(kRadix, 0);
    local_counts.combine_each([&count](const std::vector<int> &c) {
      for (int digit_idx = 0; digit_idx < kRadix; ++digit_idx) {
        count[static_cast<size_t>(digit_idx)] += c[static_cast<size_t>(digit_idx)];
      }
    });

    ComputePrefixSums(count);

    for (size_t i = n; i-- > 0;) {
      uint32_t val = temp[i];
      int digit = static_cast<int>((val >> shift) & 0xFFU);
      auto pos = static_cast<size_t>(--count[static_cast<size_t>(digit)]);
      buffer[pos] = val;
    }

    temp.swap(buffer);
  }

  tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&data, &temp](const tbb::blocked_range<size_t> &r) {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      data[i] = static_cast<int>(temp[i] ^ kSignMask);
    }
  });
}

void ChernovTRadixSortTBB::SimpleMerge(const std::vector<int> &left, const std::vector<int> &right,
                                       std::vector<int> &result) {
  result.resize(left.size() + right.size());

  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result[k++] = left[i++];
    } else {
      result[k++] = right[j++];
    }
  }
  while (i < left.size()) {
    result[k++] = left[i++];
  }
  while (j < right.size()) {
    result[k++] = right[j++];
  }
}

bool ChernovTRadixSortTBB::RunImpl() {
  auto &data = GetOutput();

  if (data.size() <= 1) {
    return true;
  }

  const size_t mid = data.size() / 2;
  std::vector<int> left(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(mid));
  std::vector<int> right(data.begin() + static_cast<std::ptrdiff_t>(mid), data.end());

  tbb::parallel_invoke([&left]() { RadixSortLSD(left); }, [&right]() { RadixSortLSD(right); });

  SimpleMerge(left, right, data);

  return true;
}

bool ChernovTRadixSortTBB::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace chernov_t_radix_sort
