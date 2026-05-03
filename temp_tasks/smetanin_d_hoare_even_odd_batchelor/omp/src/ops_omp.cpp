#include "smetanin_d_hoare_even_odd_batchelor/omp/include/ops_omp.hpp"

#include <stack>
#include <utility>
#include <vector>

#include "smetanin_d_hoare_even_odd_batchelor/common/include/common.hpp"

namespace smetanin_d_hoare_even_odd_batchelor {

namespace {

constexpr int kTaskCutoff = 1000;

int HoarePartition(std::vector<int> &arr, int lo, int hi) {
  int pivot = arr[lo + ((hi - lo) / 2)];
  int i = lo - 1;
  int j = hi + 1;
  while (true) {
    ++i;
    while (arr[i] < pivot) {
      ++i;
    }
    --j;
    while (arr[j] > pivot) {
      --j;
    }
    if (i >= j) {
      return j;
    }
    std::swap(arr[i], arr[j]);
  }
}

void OddEvenMerge(std::vector<int> &arr, int lo, int hi) {
  int n = hi - lo + 1;
  for (int step = 1; step < n; step *= 2) {
    for (int i = lo; i + step <= hi; i += step * 2) {
      if (arr[i] > arr[i + step]) {
        std::swap(arr[i], arr[i + step]);
      }
    }
  }
}

void HoarSortBatcherSeq(std::vector<int> &arr, int lo, int hi) {
  std::stack<std::pair<int, int>> stk;
  stk.emplace(lo, hi);
  while (!stk.empty()) {
    auto [l, r] = stk.top();
    stk.pop();
    if (l >= r) {
      continue;
    }
    int p = HoarePartition(arr, l, r);
    if ((p - l) > (r - p - 1)) {
      stk.emplace(l, p);
      stk.emplace(p + 1, r);
    } else {
      stk.emplace(p + 1, r);
      stk.emplace(l, p);
    }
    OddEvenMerge(arr, l, r);
  }
}

void HoarSortBatcherOMPImpl(std::vector<int> &arr, int lo, int hi) {
  if (lo >= hi) {
    return;
  }
  if (hi - lo < kTaskCutoff) {
    HoarSortBatcherSeq(arr, lo, hi);
    return;
  }
  int p = HoarePartition(arr, lo, hi);
  OddEvenMerge(arr, lo, hi);
#pragma omp task default(none) shared(arr) firstprivate(lo, p)
  HoarSortBatcherOMPImpl(arr, lo, p);
#pragma omp task default(none) shared(arr) firstprivate(hi, p)
  HoarSortBatcherOMPImpl(arr, p + 1, hi);
#pragma omp taskwait
}

}  // namespace

SmetaninDHoarSortOMP::SmetaninDHoarSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SmetaninDHoarSortOMP::ValidationImpl() {
  return true;
}

bool SmetaninDHoarSortOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SmetaninDHoarSortOMP::RunImpl() {
  auto &data = GetOutput();
  int n = static_cast<int>(data.size());
  if (n > 1) {
#pragma omp parallel default(none) shared(data, n)
#pragma omp single
    HoarSortBatcherOMPImpl(data, 0, n - 1);
  }
  return true;
}

bool SmetaninDHoarSortOMP::PostProcessingImpl() {
  return true;
}

}  // namespace smetanin_d_hoare_even_odd_batchelor
