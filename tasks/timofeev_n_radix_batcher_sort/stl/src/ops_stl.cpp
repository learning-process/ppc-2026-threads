#include "timofeev_n_radix_batcher_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

TimofeevNRadixBatcherSTL::TimofeevNRadixBatcherSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

void TimofeevNRadixBatcherSTL::CompExch(int &a, int &b, int digit) {
  int b_r = b % (digit * 10) / digit;
  int a_r = a % (digit * 10) / digit;
  if (b_r < a_r) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherSTL::BubbleSort(std::vector<int> &arr, int digit, int left, int right) {
  for (int i = left; i <= right; i++) {
    for (int j = 0; j + 1 < right - left; j++) {
      CompExch(arr[left + j], arr[left + j + 1], digit);
    }
  }
}

void TimofeevNRadixBatcherSTL::ComparR(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherSTL::OddEvenMerge(std::vector<int> &arr, int lft, int n) {
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

int TimofeevNRadixBatcherSTL::Loggo(int inputa) {
  int count = 0;
  while (inputa > 1) {
    inputa /= 2;
    count++;
  }
  return count;
}

bool TimofeevNRadixBatcherSTL::ValidationImpl() {
  return GetInput().size() >= 2;
}

bool TimofeevNRadixBatcherSTL::PreProcessingImpl() {
  return true;
}

void TimofeevNRadixBatcherSTL::PrepAux(int &n, int &m, std::vector<int> &in, int &max_x, size_t &num_threads,
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

  n_thr = std::thread::hardware_concurrency(); // from there
  num_threads = 1;
  while (num_threads * 2 <= n_thr && n / num_threads >= 4) {
    num_threads *= 2;
  }
}

void TimofeevNRadixBatcherSTL::BubbleSortAux(size_t &num_threads, int &max_x, std::vector<int> &r_in, size_t &piece,
                                             std::vector<std::thread> &threads) {
  for (size_t t_s = 0; t_s < num_threads; ++t_s) {
    threads.emplace_back([&, t_s]() {
      for (int k = 1; k <= max_x; k *= 10) {
        BubbleSort(r_in, k, static_cast<int>(piece * t_s), static_cast<int>((piece * t_s) + piece));
      }
    });
  }
}

void TimofeevNRadixBatcherSTL::OddMergeAux(std::vector<int> &r_in, size_t &piece, std::vector<std::thread> &threads,
                                           size_t &n_n) {
  for (size_t c_p = piece * 2; c_p <= n_n; c_p *= 2) {
    threads.clear();

    for (size_t i = 0; i < n_n; i += c_p) {
      threads.emplace_back([&, i]() { OddEvenMerge(r_in, static_cast<int>(i), static_cast<int>(c_p)); });
    }

    for (auto &th : threads) {
      th.join();
    }
  }
}

void TimofeevNRadixBatcherSTL::HandleSTL(size_t &num_threads, size_t &n_n, size_t &m_m, int &max_x,
                                         std::vector<int> &reff, std::vector<int> &r_in) {
  std::vector<std::thread> threads;

  size_t piece = n_n / num_threads;
  threads.reserve(num_threads);

  BubbleSortAux(num_threads, max_x, r_in, piece, threads);

  for (auto &th : threads) {
    th.join();
  }

  OddMergeAux(r_in, piece, threads, n_n);

  if (m_m != n_n) {
    r_in.resize(m_m);
  }

  threads.clear();
  size_t chunk_size = r_in.size() / num_threads;
  if (chunk_size == 0) {
    chunk_size = 1;
  }

  for (size_t t_s = 0; t_s < num_threads; ++t_s) {
    size_t start = t_s * chunk_size;
    size_t end = (t_s == num_threads - 1) ? r_in.size() : (t_s + 1) * chunk_size;

    threads.emplace_back([&, start, end]() {
      for (size_t i = start; i < end; ++i) {
        reff[i] = r_in[i];
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }
}

bool TimofeevNRadixBatcherSTL::RunImpl() {
  std::vector<int> in;
  int n = 0;
  int m = 0;
  int max_x = 0;
  size_t n_thr = 0;
  size_t num_threads = 0;

  PrepAux(n, m, in, max_x, num_threads, n_thr);

  std::vector<int> &r_in = in;
  size_t n_n = n;
  size_t m_m = m;
  std::vector<int> &reff = GetInput();

  HandleSTL(num_threads, n_n, m_m, max_x, reff, r_in);

  GetOutput() = reff;
  return true;
}

bool TimofeevNRadixBatcherSTL::PostProcessingImpl() {
  return true;
}

}  // namespace timofeev_n_radix_batcher_sort_threads
