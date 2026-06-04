#include "egashin_k_radix_simple_merge/omp/include/ops_omp.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "egashin_k_radix_simple_merge/common/include/common.hpp"
#include "egashin_k_radix_simple_merge/common/include/radix_utils.hpp"
#include "util/include/util.hpp"

namespace egashin_k_radix_simple_merge {

EgashinKRadixSimpleMergeOMP::EgashinKRadixSimpleMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool EgashinKRadixSimpleMergeOMP::ValidationImpl() {
  return true;
}

bool EgashinKRadixSimpleMergeOMP::PreProcessingImpl() {
  result_ = GetInput();
  return true;
}

bool EgashinKRadixSimpleMergeOMP::RunImpl() {
  if (result_.size() < 2) {
    return true;
  }

  int workers = radix_utils::WorkerCount(result_.size(), ppc::util::GetNumThreads());
  auto ranges = radix_utils::MakeRanges(result_.size(), workers);
  auto &result = result_;

#pragma omp parallel for default(none) shared(result, ranges, workers) num_threads(workers) schedule(static)
  for (int i = 0; i < workers; ++i) {
    radix_utils::SortRange(result, ranges[static_cast<size_t>(i)].first, ranges[static_cast<size_t>(i)].second);
  }

  auto parts = radix_utils::MakeParts(result_, ranges);
  while (parts.size() > 1) {
    size_t pair_count = parts.size() / 2;
    std::vector<std::vector<double>> next((parts.size() + 1) / 2);

#pragma omp parallel for default(none) shared(parts, next, pair_count) num_threads(workers) schedule(static)
    for (size_t i = 0; i < pair_count; ++i) {
      next[i] = radix_utils::Merge(parts[2 * i], parts[(2 * i) + 1]);
    }

    if (parts.size() % 2 != 0) {
      next.back() = std::move(parts.back());
    }
    parts = std::move(next);
  }

  result_ = std::move(parts.front());
  return true;
}

bool EgashinKRadixSimpleMergeOMP::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace egashin_k_radix_simple_merge
