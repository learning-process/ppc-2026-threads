#include "egashin_k_radix_simple_merge/stl/include/ops_stl.hpp"

#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "egashin_k_radix_simple_merge/common/include/common.hpp"
#include "egashin_k_radix_simple_merge/common/include/radix_utils.hpp"
#include "util/include/util.hpp"

namespace egashin_k_radix_simple_merge {

EgashinKRadixSimpleMergeSTL::EgashinKRadixSimpleMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool EgashinKRadixSimpleMergeSTL::ValidationImpl() {
  return true;
}

bool EgashinKRadixSimpleMergeSTL::PreProcessingImpl() {
  result_ = GetInput();
  return true;
}

bool EgashinKRadixSimpleMergeSTL::RunImpl() {
  if (result_.size() < 2) {
    return true;
  }

  const int workers = radix_utils::WorkerCount(result_.size(), ppc::util::GetNumThreads());
  const auto ranges = radix_utils::MakeRanges(result_.size(), workers);

  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(workers));
  for (int i = 0; i < workers; ++i) {
    threads.emplace_back([&, index = static_cast<size_t>(i)] {
      radix_utils::SortRange(result_, ranges[index].first, ranges[index].second);
    });
  }
  for (auto &thread : threads) {
    thread.join();
  }

  auto parts = radix_utils::MakeParts(result_, ranges);
  while (parts.size() > 1) {
    const size_t pair_count = parts.size() / 2;
    std::vector<std::vector<double>> next((parts.size() + 1) / 2);

    const int merge_workers = radix_utils::WorkerCount(pair_count, workers);
    const auto merge_ranges = radix_utils::MakeRanges(pair_count, merge_workers);

    threads.clear();
    threads.reserve(static_cast<size_t>(merge_workers));
    for (int i = 0; i < merge_workers; ++i) {
      threads.emplace_back([&, index = static_cast<size_t>(i)] {
        for (size_t j = merge_ranges[index].first; j < merge_ranges[index].second; ++j) {
          next[j] = radix_utils::Merge(parts[2 * j], parts[(2 * j) + 1]);
        }
      });
    }
    for (auto &thread : threads) {
      thread.join();
    }

    if (parts.size() % 2 != 0) {
      next.back() = std::move(parts.back());
    }
    parts = std::move(next);
  }

  result_ = std::move(parts.front());
  return true;
}

bool EgashinKRadixSimpleMergeSTL::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace egashin_k_radix_simple_merge
