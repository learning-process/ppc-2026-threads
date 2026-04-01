#include "timofeev_n_radix_batcher_sort/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <utility>
#include <vector>

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

TimofeevNRadixBatcherTBB::TimofeevNRadixBatcherTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

void TimofeevNRadixBatcherTBB::CompExch(int &a, int &b, int digit) {
  int b_r = b % (digit * 10) / digit;
  int a_r = a % (digit * 10) / digit;
  if (b_r < a_r) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherTBB::BubbleSort(std::vector<int> &arr, int digit, int left, int right) {
  for (int i = left; i <= right; i++) {
    for (int j = 0; j + 1 < right - left; j++) {
      CompExch(arr[left + j], arr[left + j + 1], digit);
    }
  }
}

void TimofeevNRadixBatcherTBB::ComparR(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherTBB::OddEvenMerge(std::vector<int> &arr, int lft, int n) {
  if (n <= 1) {
    return;
  }

  int otstup = n / 2;
  for (int i = 0; i < otstup; i += 1) {
    if (arr[lft + i] > arr[lft + otstup + i]) {
      std::swap(arr[lft + i], arr[lft + otstup + i]);
    }
  }

  for (otstup = n / 4; otstup > 0; otstup /= 2) {
    int h = otstup * 2;
    for (int start = otstup; start + otstup < n; start += h) {
      for (int i = 0; i < otstup; i += 1) {
        ComparR(arr[lft + start + i], arr[lft + start + i + otstup]);
      }
    }
  }
}

int TimofeevNRadixBatcherTBB::Loggo(int inputa) {
  int count = 0;
  while (inputa > 1) {
    inputa /= 2;
    count++;
  }
  return count;
}

bool TimofeevNRadixBatcherTBB::ValidationImpl() {
  return GetInput().size() >= 2;
}

bool TimofeevNRadixBatcherTBB::PreProcessingImpl() {
  return true;
}

void TimofeevNRadixBatcherTBB::PrepAux(int &n, int &m, std::vector<int> &in, int &max_x, size_t &num_threads,
                                       size_t &n_thr) {
  in = GetInput();
  n = static_cast<int>(in.size());
  m = n;

  while (n % 2 == 0) {
    n /= 2;
  }
  if (n > 1) {
    n = static_cast<int>(in.size());
    int p = 1;
    while (p < n) {
      p *= 2;
    }
    n = p;
  } else {
    n = m;
  }

  max_x = *(std::ranges::max_element(in.begin(), in.end()));
  if (n != m) {
    in.resize(n, max_x);
  }
  // tbb
  n_thr = tbb::this_task_arena::max_concurrency();
  num_threads = 1;
  while (num_threads * 2 <= n_thr && n / num_threads >= 4) {
    num_threads *= 2;
  }
}

void TimofeevNRadixBatcherTBB::HandleTBB(size_t &num_threads, size_t &n_n, size_t &m_m, int &max_x, std::vector<int> &reff, std::vector<int> &r_in) {
  tbb::task_arena arena(static_cast<int>(num_threads));
  arena.execute([&] {
    size_t piece = 0;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_threads), [&](const tbb::blocked_range<size_t> &range) {
      size_t t_n = range.begin();
      piece = n_n / num_threads;
      for (int k = 1; k <= max_x; k *= 10) {
        BubbleSort(r_in, k, static_cast<int>(piece * t_n), static_cast<int>((piece * t_n) + piece));
      }
    });

    size_t c_p = piece * 2;
    for (; c_p <= n_n; c_p *= 2) {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, n_n, c_p), [&](const tbb::blocked_range<size_t> &range) {
        for (size_t i = range.begin(); i < range.end(); i += c_p) {
          OddEvenMerge(r_in, static_cast<int>(i), static_cast<int>(c_p));
        }
      });
    }

    if (m_m != n_n && tbb::this_task_arena::current_thread_index() == 0) {
      r_in.resize(m_m);
    }

    tbb::parallel_for(tbb::blocked_range<size_t>(0, r_in.size()), [&](const tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i < range.end(); ++i) {
        reff[i] = r_in[i];
      }
    });
  });
}

bool TimofeevNRadixBatcherTBB::RunImpl() {
  std::vector<int> in;
  int n;
  int m;
  int max_x;
  size_t n_thr;  // = tbb::this_task_arena::max_concurrency();
  size_t num_threads;

  PrepAux(n, m, in, max_x, num_threads, n_thr);

  std::vector<int> &r_in = in;
  size_t n_n = n;
  size_t m_m = m;
  std::vector<int> &reff = GetInput();

  HandleTBB(num_threads, n_n, m_m, max_x, reff, r_in);
  // tbb
  GetOutput() = reff;
  return true;
}

bool TimofeevNRadixBatcherTBB::PostProcessingImpl() {
  return true;
}

}  // namespace timofeev_n_radix_batcher_sort_threads
