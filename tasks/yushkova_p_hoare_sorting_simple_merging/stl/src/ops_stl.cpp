#include "yushkova_p_hoare_sorting_simple_merging/stl/include/ops_stl.hpp"

#include <algorithm>
#include <exception>
#include <iterator>
#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

YushkovaPHoareSortingSimpleMergingSTL::YushkovaPHoareSortingSimpleMergingSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

int YushkovaPHoareSortingSimpleMergingSTL::HoarePartition(std::vector<int> &values, int left, int right) {
  const int pivot = values[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    ++i;
    while (values[i] < pivot) {
      ++i;
    }

    --j;
    while (values[j] > pivot) {
      --j;
    }

    if (i >= j) {
      return j;
    }

    std::swap(values[i], values[j]);
  }
}

void YushkovaPHoareSortingSimpleMergingSTL::HoareQuickSort(std::vector<int> &values, int left, int right) {
  std::stack<std::pair<int, int>> ranges;
  ranges.emplace(left, right);

  while (!ranges.empty()) {
    auto [current_left, current_right] = ranges.top();
    ranges.pop();

    if (current_left >= current_right) {
      continue;
    }

    const int partition_index = HoarePartition(values, current_left, current_right);

    if ((partition_index - current_left) > (current_right - (partition_index + 1))) {
      ranges.emplace(current_left, partition_index);
      ranges.emplace(partition_index + 1, current_right);
    } else {
      ranges.emplace(partition_index + 1, current_right);
      ranges.emplace(current_left, partition_index);
    }
  }
}

std::vector<int> YushkovaPHoareSortingSimpleMergingSTL::SimpleMerge(const std::vector<int> &left,
                                                                    const std::vector<int> &right) {
  std::vector<int> merged;
  merged.reserve(left.size() + right.size());
  std::ranges::merge(left, right, std::back_inserter(merged));
  return merged;
}

void YushkovaPHoareSortingSimpleMergingSTL::SortHalfIfNeeded(std::vector<int> &values) {
  if (values.size() > 1) {
    HoareQuickSort(values, 0, static_cast<int>(values.size()) - 1);
  }
}

void YushkovaPHoareSortingSimpleMergingSTL::SortHalvesSequential(std::vector<int> &left, std::vector<int> &right) {
  SortHalfIfNeeded(left);
  SortHalfIfNeeded(right);
}

void YushkovaPHoareSortingSimpleMergingSTL::SortHalvesParallel(std::vector<int> &left, std::vector<int> &right) {
  std::exception_ptr exception_ptr;
  std::thread left_worker([&] {
    try {
      SortHalfIfNeeded(left);
    } catch (...) {
      exception_ptr = std::current_exception();
    }
  });

  try {
    SortHalfIfNeeded(right);
  } catch (...) {
    if (!exception_ptr) {
      exception_ptr = std::current_exception();
    }
  }

  left_worker.join();

  if (exception_ptr) {
    std::rethrow_exception(exception_ptr);
  }
}

bool YushkovaPHoareSortingSimpleMergingSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool YushkovaPHoareSortingSimpleMergingSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool YushkovaPHoareSortingSimpleMergingSTL::RunImpl() {
  std::vector<int> &values = GetOutput();
  if (values.size() <= 1) {
    return true;
  }

  const auto middle = static_cast<std::vector<int>::difference_type>(values.size() / 2);
  std::vector<int> left(values.begin(), values.begin() + middle);
  std::vector<int> right(values.begin() + middle, values.end());

  const int concurrency = std::max(1, ppc::util::GetNumThreads());
  if (concurrency == 1) {
    SortHalvesSequential(left, right);
  } else {
    SortHalvesParallel(left, right);
  }

  values = SimpleMerge(left, right);
  return std::ranges::is_sorted(values);
}

bool YushkovaPHoareSortingSimpleMergingSTL::PostProcessingImpl() {
  return !GetOutput().empty() && std::ranges::is_sorted(GetOutput());
}

}  // namespace yushkova_p_hoare_sorting_simple_merging
