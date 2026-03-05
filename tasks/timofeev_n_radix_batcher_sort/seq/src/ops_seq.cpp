#include "timofeev_n_radix_batcher_sort/seq/include/ops_seq.hpp"

#include <climits>
#include <numeric>

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

TimofeevNRadixBatcherSEQ::TimofeevNRadixBatcherSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

void TimofeevNRadixBatcherSEQ::CompExch(int &a, int &b, int digit) {
  int bR = b % (digit * 10) / digit;
  int aR = a % (digit * 10) / digit;
  if (bR < aR) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherSEQ::BubbleSort(std::vector<int> &arr, int digit, int left, int right) {
  for (int i = left; i <= right; i++) {
    for (int j = 0; j + 1 < right - left; j++) {
      CompExch(arr[left + j], arr[left + j + 1], digit);
    }
  }
}

void TimofeevNRadixBatcherSEQ::ComparR(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherSEQ::OddEvenMerge(std::vector<int> &arr, size_t lft, size_t n) {
  if (n <= 1) {
    return;
  }

  size_t otstup = n / 2;
  for (size_t i = 0; i < otstup; i += 1) {
    if (arr[lft + i] > arr[lft + otstup + i]) {
      std::swap(arr[lft + i], arr[lft + otstup + i]);
    }
  }

  for (otstup = n / 4; otstup > 0; otstup /= 2) {
    size_t h = otstup * 2;
    for (size_t start = otstup; start + otstup < n; start += h) {
      for (size_t i = 0; i < otstup; i += 1) {
        ComparR(arr[lft + start + i], arr[lft + start + i + otstup]);
      }
    }
  }
}

int TimofeevNRadixBatcherSEQ::Loggo(int inputa) {
  int count = 0;
  while (inputa > 1) {
    inputa /= 2;
    count++;
  }
  return count;
}

bool TimofeevNRadixBatcherSEQ::ValidationImpl() {
  if (GetInput().size() < 2) {
    return false;
  }
  return true;
}

bool TimofeevNRadixBatcherSEQ::PreProcessingImpl() {
  return true;
}

bool TimofeevNRadixBatcherSEQ::RunImpl() {
  std::vector<int> in = GetInput();
  int n = in.size();
  int m = n;
  while (n % 2 == 0) {
    n /= 2;
  }
  if (n > 1) {
    n = in.size();
    int p = 1;
    while (p < n) {
      p *= 2;
    }
    n = p;
  } else {
    n = m;
  }
  int maxX = *(std::max_element(in.begin(), in.end()));
  if (n != m) {
    in.resize(n, maxX);
  }
  // std::cout << "\n\n" << n << " " << in.size() << " " << "1\n\n";
  for (int i = 0; i < static_cast<int>(in.size()); i += static_cast<int>(in.size() / 2)) {
    for (int k = 1; k <= maxX; k *= 10) {
      BubbleSort(in, k, i, i + static_cast<int>(in.size() / 2));
    }
  }
  OddEvenMerge(in, 0, in.size());
  // std::cout << "\n\n" << n << " " << in.size() << " " << "2\n\n";
  if (m != n) {
    in.resize(m);
  }
  for (int i = 0; i < static_cast<int>(in.size()); i++) {
    GetInput()[i] = in[i];
  }
  GetOutput() = GetInput();
  return true;
}

bool TimofeevNRadixBatcherSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace timofeev_n_radix_batcher_sort_threads
