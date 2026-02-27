#include "vasiliev_m_shell_sort_batcher_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

VasilievMShellSortBatcherMergeSEQ::VasilievMShellSortBatcherMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool VasilievMShellSortBatcherMergeSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool VasilievMShellSortBatcherMergeSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool VasilievMShellSortBatcherMergeSEQ::RunImpl() {
  auto &vec = GetInput();
  const size_t n = vec.size();

  if (vec.empty()) {
    return false;
  }

  size_t mid = n / 2;

  std::vector<ValType> l(vec.begin(), vec.begin() + static_cast<std::ptrdiff_t>(mid));
  std::vector<ValType> r(vec.begin() + static_cast<std::ptrdiff_t>(mid), vec.end());

  ShellSort(l);
  ShellSort(r);

  std::vector<ValType> batcher_merged_vec = BatcherMerge(l, r);

  GetOutput() = batcher_merged_vec;
  return true;
}

void VasilievMShellSortBatcherMergeSEQ::ShellSort(std::vector<ValType> &vec) {
  size_t n = vec.size();
  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; i++) {
      ValType tmp = vec[i];
      size_t j = i;
      while (j >= gap && vec[j - gap] > tmp) {
        vec[j] = vec[j - gap];
        j -= gap;
      }
      vec[j] = tmp;
    }
  }
}

std::vector<ValType> VasilievMShellSortBatcherMergeSEQ::BatcherMerge(std::vector<ValType> &l, std::vector<ValType> &r) {
  std::vector<ValType> even_l;
  std::vector<ValType> odd_l;
  std::vector<ValType> even_r;
  std::vector<ValType> odd_r;

  SplitEvenOdd(l, even_l, odd_l);
  SplitEvenOdd(r, even_r, odd_r);

  std::vector<ValType> even = Merge(even_l, even_r);
  std::vector<ValType> odd = Merge(odd_l, odd_r);

  std::vector<ValType> res;
  res.reserve(l.size() + r.size());

  for (size_t i = 0; i < even.size() || i < odd.size(); i++) {
    if (i < even.size()) {
      res.push_back(even[i]);
    }
    if (i < odd.size()) {
      res.push_back(odd[i]);
    }
  }

  for (size_t i = 1; i + 1 < res.size(); i += 2) {
    if (res[i] > res[i + 1]) {
      std::swap(res[i], res[i + 1]);
    }
  }

  return res;
}

void VasilievMShellSortBatcherMergeSEQ::SplitEvenOdd(std::vector<ValType> &vec, std::vector<ValType> &even,
                                                     std::vector<ValType> &odd) {
  even.reserve(even.size() + (vec.size() / 2) + 1);
  odd.reserve(odd.size() + (vec.size() / 2));

  for (size_t i = 0; i < vec.size(); i += 2) {
    even.push_back(vec[i]);
    if (i + 1 < vec.size()) {
      odd.push_back(vec[i + 1]);
    }
  }
}

std::vector<ValType> VasilievMShellSortBatcherMergeSEQ::Merge(std::vector<ValType> &a, std::vector<ValType> &b) {
  std::vector<ValType> merged;
  size_t i = 0;
  size_t j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      merged.push_back(a[i++]);
    } else {
      merged.push_back(b[j++]);
    }
  }

  while (i < a.size()) {
    merged.push_back(a[i++]);
  }
  while (j < b.size()) {
    merged.push_back(b[j++]);
  }

  return merged;
}

bool VasilievMShellSortBatcherMergeSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace vasiliev_m_shell_sort_batcher_merge
