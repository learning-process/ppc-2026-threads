#include "shekhirev_v_hoare_batcher_sort_seq/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <climits>
#include <utility>
#include <vector>

namespace shekhirev_v_hoare_batcher_sort_seq {

namespace {
void HoareSort(std::vector<int> &arr, int low, int high) {
  if (low >= high) {
    return;
  }
  int pivot = arr[low + (high - low) / 2];
  int i = low - 1;
  int j = high + 1;
  while (true) {
    do {
      i++;
    } while (arr[i] < pivot);
    do {
      j--;
    } while (arr[j] > pivot);
    if (i >= j) {
      break;
    }
    std::swap(arr[i], arr[j]);
  }
  HoareSort(arr, low, j);
  HoareSort(arr, j + 1, high);
}
}  // namespace

ShekhirevHoareBatcherSortOMP::ShekhirevHoareBatcherSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ShekhirevHoareBatcherSortOMP::ValidationImpl() {
  return true;
}

bool ShekhirevHoareBatcherSortOMP::PreProcessingImpl() {
  input_ = GetInput();
  return true;
}

bool ShekhirevHoareBatcherSortOMP::RunImpl() {
  int n = static_cast<int>(input_.size());
  if (n <= 1) {
    res_ = input_;
    return true;
  }

  int n_pow2 = 1;
  while (n_pow2 < n) {
    n_pow2 *= 2;
  }

  std::vector<int> a(n_pow2, INT_MAX);
  for (int i = 0; i < n; i++) {
    a[i] = input_[i];
  }

  int num_threads = omp_get_max_threads();
  if (num_threads <= 0) {
    num_threads = 1;
  }

  int p_threads = 1;
  while (p_threads * 2 <= num_threads && p_threads * 2 <= n_pow2) {
    p_threads *= 2;
  }

  int chunk_size = n_pow2 / p_threads;

#pragma omp parallel for
  for (int i = 0; i < p_threads; i++) {
    HoareSort(a, i * chunk_size, (i + 1) * chunk_size - 1);
  }

  for (int p = chunk_size; p < n_pow2; p *= 2) {
    for (int k = p; k >= 1; k /= 2) {
      int start = (k == p) ? 0 : k;
      int num_blocks = (start == 0) ? (n_pow2 / (2 * k)) : (n_pow2 / (2 * k) - 1);
      int total_pairs = num_blocks * k;

#pragma omp parallel for
      for (int step = 0; step < total_pairs; ++step) {
        int i = step % k;
        int b = step / k;
        int j = start + b * (k * 2);

        int idx1 = i + j;
        int idx2 = i + j + k;

        if (idx1 / (p * 2) == idx2 / (p * 2)) {
          if (a[idx1] > a[idx2]) {
            std::swap(a[idx1], a[idx2]);
          }
        }
      }
    }
  }

  res_.assign(a.begin(), a.begin() + n);
  return true;
}

bool ShekhirevHoareBatcherSortOMP::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace shekhirev_v_hoare_batcher_sort_seq
