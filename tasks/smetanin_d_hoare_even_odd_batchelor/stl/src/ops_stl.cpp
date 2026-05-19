#include "smetanin_d_hoare_even_odd_batchelor/stl/include/ops_stl.hpp"

#include <algorithm>
#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "smetanin_d_hoare_even_odd_batchelor/common/include/common.hpp"
#include "util/include/util.hpp"

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

void HoarSortBatcherSTLImpl(std::vector<int> &arr, int lo, int hi, int num_threads) {
  if (lo >= hi) {
    return;
  }
  if (hi - lo < kTaskCutoff || num_threads <= 1) {
    HoarSortBatcherSeq(arr, lo, hi);
    return;
  }
  const int p = HoarePartition(arr, lo, hi);
  OddEvenMerge(arr, lo, hi);
  std::thread right([&arr, p, hi]() { HoarSortBatcherSeq(arr, p + 1, hi); });
  HoarSortBatcherSeq(arr, lo, p);
  right.join();
}

}  // namespace

SmetaninDHoarSortSTL::SmetaninDHoarSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SmetaninDHoarSortSTL::ValidationImpl() {
  return true;
}

bool SmetaninDHoarSortSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SmetaninDHoarSortSTL::RunImpl() {
  auto &data = GetOutput();
  int n = static_cast<int>(data.size());
  if (n > 1) {
    const int num_threads = std::max(1, ppc::util::GetNumThreads());
    HoarSortBatcherSTLImpl(data, 0, n - 1, num_threads);
  }
  return true;
}

bool SmetaninDHoarSortSTL::PostProcessingImpl() {
  return true;
}

}  // namespace smetanin_d_hoare_even_odd_batchelor
