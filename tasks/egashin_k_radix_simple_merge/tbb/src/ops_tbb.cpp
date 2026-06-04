#include "egashin_k_radix_simple_merge/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/parallel_for.h>

#include <cstddef>
#include <utility>
#include <vector>

#include "egashin_k_radix_simple_merge/common/include/common.hpp"
#include "egashin_k_radix_simple_merge/common/include/radix_utils.hpp"
#include "util/include/util.hpp"

namespace egashin_k_radix_simple_merge {

EgashinKRadixSimpleMergeTBB::EgashinKRadixSimpleMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool EgashinKRadixSimpleMergeTBB::ValidationImpl() {
  return true;
}

bool EgashinKRadixSimpleMergeTBB::PreProcessingImpl() {
  result_ = GetInput();
  return true;
}

bool EgashinKRadixSimpleMergeTBB::RunImpl() {
  if (result_.size() < 2) {
    return true;
  }

  const int workers = radix_utils::WorkerCount(result_.size(), ppc::util::GetNumThreads());
  oneapi::tbb::global_control control(oneapi::tbb::global_control::max_allowed_parallelism, workers);

  const auto ranges = radix_utils::MakeRanges(result_.size(), workers);
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, ranges.size()),
                            [&](const oneapi::tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      radix_utils::SortRange(result_, ranges[i].first, ranges[i].second);
    }
  });

  auto parts = radix_utils::MakeParts(result_, ranges);
  while (parts.size() > 1) {
    const size_t pair_count = parts.size() / 2;
    std::vector<std::vector<double>> next((parts.size() + 1) / 2);

    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, pair_count),
                              [&](const oneapi::tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i != range.end(); ++i) {
        next[i] = radix_utils::Merge(parts[2 * i], parts[(2 * i) + 1]);
      }
    });

    if (parts.size() % 2 != 0) {
      next.back() = std::move(parts.back());
    }
    parts = std::move(next);
  }

  result_ = std::move(parts.front());
  return true;
}

bool EgashinKRadixSimpleMergeTBB::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace egashin_k_radix_simple_merge
