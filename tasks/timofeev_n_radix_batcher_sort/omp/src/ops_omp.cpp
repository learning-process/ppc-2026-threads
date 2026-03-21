#include "timofeev_n_radix_batcher_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <utility>
#include <vector>

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

TimofeevNRadixBatcherOMP::TimofeevNRadixBatcherOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

void TimofeevNRadixBatcherOMP::CompExch(int &a, int &b, int digit) {
  int b_r = b % (digit * 10) / digit;
  int a_r = a % (digit * 10) / digit;
  if (b_r < a_r) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherOMP::BubbleSort(std::vector<int> &arr, int digit, int left, int right) {
  for (int i = left; i <= right; i++) {
    for (int j = 0; j + 1 < right - left; j++) {
      CompExch(arr[left + j], arr[left + j + 1], digit);
    }
  }
}

void TimofeevNRadixBatcherOMP::ComparR(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherOMP::OddEvenMerge(std::vector<int> &arr, int lft, int n) {
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

int TimofeevNRadixBatcherOMP::Loggo(int inputa) {
  int count = 0;
  while (inputa > 1) {
    inputa /= 2;
    count++;
  }
  return count;
}

bool TimofeevNRadixBatcherOMP::ValidationImpl() {
  return GetInput().size() >= 2;
}

bool TimofeevNRadixBatcherOMP::PreProcessingImpl() {
  return true;
}

bool TimofeevNRadixBatcherOMP::RunImpl() {
  std::vector<int> in = GetInput();
  int n = static_cast<int>(in.size());
  int m = n;
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
  int max_x = *(std::ranges::max_element(in.begin(), in.end()));
  if (n != m) {
    in.resize(n, max_x);
  }
  size_t n_thr = omp_get_max_threads();
  size_t num_threads = 1;
  while (num_threads * 2 <= n_thr && n / num_threads >= 4) {
    num_threads *= 2;
  }
  std::vector<int> &r_in = in;
  size_t n_n = n;
  size_t m_m = m;
  std::vector<int> &reff = GetInput();
#pragma omp parallel num_threads(num_threads) default(none) shared(r_in, max_x, n_n, num_threads, m_m, reff)
  {
    size_t t_n = omp_get_thread_num();
    size_t piece = n_n / num_threads;
    for (int k = 1; k <= max_x; k *= 10) {
      BubbleSort(r_in, k, static_cast<int>(piece * t_n), static_cast<int>((piece * t_n) + piece));  // [left; right)
    }
#pragma omp barrier
    size_t c_p = piece * 2;
    for (; c_p <= n_n; c_p *= 2) {
#pragma omp for
      for (size_t i = 0; i < n_n; i += c_p) {
        OddEvenMerge(r_in, static_cast<int>(i), static_cast<int>(c_p));
      }
#pragma omp barrier
    }
    if (t_n == 0 && m_m != n_n) {
      r_in.resize(m_m);
    }
#pragma omp barrier
#pragma omp for
    for (size_t i = 0; i < r_in.size(); i++) {
      reff[i] = r_in[i];
    }
#pragma omp barrier
  }

  GetOutput() = reff;
  return true;
}

bool TimofeevNRadixBatcherOMP::PostProcessingImpl() {
  return true;
}

}  // namespace timofeev_n_radix_batcher_sort_threads
