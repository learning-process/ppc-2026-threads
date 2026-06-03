#include "shelenkova_m_shell_sort_simple_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <thread>
#include <vector>

#include "shelenkova_m_shell_sort_simple_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shelenkova_m_shell_sort_simple_merge {

namespace {

void ShellSort(std::vector<int>::iterator begin, std::vector<int>::iterator end) {
  const size_t n = static_cast<size_t>(std::distance(begin, end));
  if (n <= 1) {
    return;
  }

  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; ++i) {
      int temp = *(begin + static_cast<std::ptrdiff_t>(i));
      size_t j = i;
      while (j >= gap && *(begin + static_cast<std::ptrdiff_t>(j - gap)) > temp) {
        *(begin + static_cast<std::ptrdiff_t>(j)) = *(begin + static_cast<std::ptrdiff_t>(j - gap));
        j -= gap;
      }
      *(begin + static_cast<std::ptrdiff_t>(j)) = temp;
    }
  }
}

}  // namespace

ShelenkovaMShellSortSimpleMergeSTL::ShelenkovaMShellSortSimpleMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool ShelenkovaMShellSortSimpleMergeSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool ShelenkovaMShellSortSimpleMergeSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool ShelenkovaMShellSortSimpleMergeSTL::RunImpl() {
  std::vector<int> &data = GetOutput();
  const size_t n = data.size();
  if (n <= 1) {
    return true;
  }

  const int num_threads = std::max(1, ppc::util::GetNumThreads());
  const size_t thread_count = std::min<size_t>(static_cast<size_t>(num_threads), n);
  const size_t chunk = n / thread_count;
  const size_t rem = n % thread_count;

  std::vector<size_t> offsets(thread_count + 1);
  size_t curr = 0;
  for (size_t i = 0; i < thread_count; ++i) {
    offsets[i] = curr;
    curr += chunk + (i < rem ? 1 : 0);
  }
  offsets[thread_count] = n;

  {
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    for (size_t i = 0; i < thread_count; ++i) {
      threads.emplace_back([&data, &offsets, i]() {
        ShellSort(data.begin() + static_cast<std::ptrdiff_t>(offsets[i]),
                  data.begin() + static_cast<std::ptrdiff_t>(offsets[i + 1]));
      });
    }
    for (auto &t : threads) {
      t.join();
    }
  }

  for (size_t step = 1; step < thread_count; ++step) {
    std::inplace_merge(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(offsets[step]),
                       data.begin() + static_cast<std::ptrdiff_t>(offsets[step + 1]));
  }

  return std::ranges::is_sorted(data);
}

bool ShelenkovaMShellSortSimpleMergeSTL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace shelenkova_m_shell_sort_simple_merge
