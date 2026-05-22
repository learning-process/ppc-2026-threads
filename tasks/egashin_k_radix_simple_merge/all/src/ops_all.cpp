#include "egashin_k_radix_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>

#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "egashin_k_radix_simple_merge/common/include/common.hpp"
#include "egashin_k_radix_simple_merge/common/include/radix_utils.hpp"
#include "util/include/util.hpp"

namespace egashin_k_radix_simple_merge {

namespace {

std::vector<int> MakeCounts(size_t size, int rank_count) {
  std::vector<int> counts(static_cast<size_t>(rank_count), 0);
  const size_t base = size / static_cast<size_t>(rank_count);
  const size_t extra = size % static_cast<size_t>(rank_count);

  for (int rank = 0; rank < rank_count; ++rank) {
    const auto rank_index = static_cast<size_t>(rank);
    counts[rank_index] = static_cast<int>(base + (std::cmp_less(rank, extra) ? size_t{1} : size_t{0}));
  }
  return counts;
}

std::vector<int> MakeDispls(const std::vector<int> &counts) {
  std::vector<int> displs(counts.size(), 0);
  for (size_t i = 1; i < counts.size(); ++i) {
    displs[i] = displs[i - 1] + counts[i - 1];
  }
  return displs;
}

void SortLocal(std::vector<double> &data) {
  if (data.size() < 2) {
    return;
  }

  int workers = radix_utils::WorkerCount(data.size(), ppc::util::GetNumThreads());
  auto ranges = radix_utils::MakeRanges(data.size(), workers);

#pragma omp parallel for default(none) shared(data, ranges, workers) num_threads(workers) schedule(static)
  for (int i = 0; i < workers; ++i) {
    radix_utils::SortRange(data, ranges[static_cast<size_t>(i)].first, ranges[static_cast<size_t>(i)].second);
  }

  auto parts = radix_utils::MakeParts(data, ranges);
  while (parts.size() > 1) {
    size_t pair_count = parts.size() / 2;
    std::vector<std::vector<double>> next((parts.size() + 1) / 2);

#pragma omp parallel for default(none) shared(parts, next, pair_count, workers) num_threads(workers) schedule(static)
    for (size_t i = 0; i < pair_count; ++i) {
      next[i] = radix_utils::Merge(parts[2 * i], parts[(2 * i) + 1]);
    }

    if (parts.size() % 2 != 0) {
      next.back() = std::move(parts.back());
    }
    parts = std::move(next);
  }

  data = std::move(parts.front());
}

std::vector<double> MergeGathered(const std::vector<double> &data, const std::vector<int> &counts,
                                  const std::vector<int> &displs) {
  std::vector<std::vector<double>> parts(counts.size());
  for (size_t rank = 0; rank < counts.size(); ++rank) {
    const auto begin = data.begin() + displs[rank];
    parts[rank] = std::vector<double>(begin, begin + counts[rank]);
  }

  while (parts.size() > 1) {
    const size_t pair_count = parts.size() / 2;
    std::vector<std::vector<double>> next((parts.size() + 1) / 2);

    for (size_t i = 0; i < pair_count; ++i) {
      next[i] = radix_utils::Merge(parts[2 * i], parts[(2 * i) + 1]);
    }

    if (parts.size() % 2 != 0) {
      next.back() = std::move(parts.back());
    }
    parts = std::move(next);
  }

  return parts.empty() ? std::vector<double>{} : std::move(parts.front());
}

}  // namespace

EgashinKRadixSimpleMergeALL::EgashinKRadixSimpleMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool EgashinKRadixSimpleMergeALL::ValidationImpl() {
  return GetInput().size() <= static_cast<size_t>(std::numeric_limits<int>::max());
}

bool EgashinKRadixSimpleMergeALL::PreProcessingImpl() {
  result_ = GetInput();
  return true;
}

bool EgashinKRadixSimpleMergeALL::RunImpl() {
  int rank = 0;
  int rank_count = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &rank_count);

  const auto counts = MakeCounts(result_.size(), rank_count);
  const auto displs = MakeDispls(counts);
  const int local_count = counts[static_cast<size_t>(rank)];
  const int local_shift = displs[static_cast<size_t>(rank)];

  std::vector<double> local(result_.begin() + local_shift, result_.begin() + local_shift + local_count);
  SortLocal(local);

  std::vector<double> gathered;
  if (rank == 0) {
    gathered = result_;
  }

  MPI_Gatherv(local.data(), local_count, MPI_DOUBLE, gathered.data(), counts.data(), displs.data(), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  if (rank == 0) {
    result_ = MergeGathered(gathered, counts, displs);
  } else {
    result_.resize(GetInput().size());
  }

  MPI_Bcast(result_.data(), static_cast<int>(result_.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return true;
}

bool EgashinKRadixSimpleMergeALL::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace egashin_k_radix_simple_merge
