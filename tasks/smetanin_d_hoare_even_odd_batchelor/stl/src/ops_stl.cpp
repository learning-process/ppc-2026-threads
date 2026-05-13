#include "smetanin_d_hoare_even_odd_batchelor/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <mutex>
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

  std::vector<std::pair<int, int>> cur;
  cur.emplace_back(lo, hi);

  while (!cur.empty()) {
    std::vector<std::pair<int, int>> next;
    std::mutex next_mtx;
    const std::size_t cur_sz = cur.size();

    auto process_range = [&](std::size_t idx) {
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
      std::lock_guard<std::mutex> lk(next_mtx);
      next.push_back({l, p});
      next.push_back({p + 1, r});
    };

    const int workers = std::clamp(num_threads, 1, std::max<int>(1, static_cast<int>(cur_sz)));

    if (workers <= 1 || cur_sz <= 1) {
      for (std::size_t i = 0; i < cur_sz; ++i) {
        process_range(i);
      }
    } else {
      std::vector<std::thread> threads;
      std::atomic<std::size_t> cursor{0};
      threads.reserve(static_cast<std::size_t>(workers));
      for (int w = 0; w < workers; ++w) {
        threads.emplace_back([&]() {
          while (true) {
            const std::size_t i = cursor.fetch_add(1);
            if (i >= cur_sz) {
              break;
            }
            process_range(i);
          }
        });
      }
      for (auto &t : threads) {
        t.join();
      }
    }

    cur = std::move(next);
  }
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
