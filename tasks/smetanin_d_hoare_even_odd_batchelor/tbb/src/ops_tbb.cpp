#include "smetanin_d_hoare_even_odd_batchelor/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/concurrent_vector.h>
#include <oneapi/tbb/parallel_for.h>

#include <cstddef>
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

void HoarSortBatcherTBBImpl(std::vector<int> &arr, int lo, int hi) {
  if (lo >= hi) {
    return;
  }

  std::vector<std::pair<int, int>> cur;
  cur.emplace_back(lo, hi);

  while (!cur.empty()) {
    tbb::concurrent_vector<std::pair<int, int>> next;
    const std::size_t cur_sz = cur.size();
    tbb::parallel_for(static_cast<std::size_t>(0), cur_sz, [&](std::size_t idx) {
      const int l = cur[idx].first;
      const int r = cur[idx].second;
      if (l >= r) {
        return;
      }
      if (r - l < kTaskCutoff) {
        HoarSortBatcherSeq(arr, l, r);
        return;
      }
      const int p = HoarePartition(arr, l, r);
      OddEvenMerge(arr, l, r);
      next.push_back({l, p});
      next.push_back({p + 1, r});
    });
    cur.assign(next.begin(), next.end());
  }
}

}  // namespace

SmetaninDHoarSortTBB::SmetaninDHoarSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SmetaninDHoarSortTBB::ValidationImpl() {
  return true;
}

bool SmetaninDHoarSortTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SmetaninDHoarSortTBB::RunImpl() {
  auto &data = GetOutput();
  int n = static_cast<int>(data.size());
  if (n > 1) {
    HoarSortBatcherTBBImpl(data, 0, n - 1);
  }
  return true;
}

bool SmetaninDHoarSortTBB::PostProcessingImpl() {
  return true;
}

}  // namespace smetanin_d_hoare_even_odd_batchelor
