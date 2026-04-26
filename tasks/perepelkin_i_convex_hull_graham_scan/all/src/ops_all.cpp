#include "perepelkin_i_convex_hull_graham_scan/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "util/include/util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

PerepelkinIConvexHullGrahamScanALL::PerepelkinIConvexHullGrahamScanALL(const InType &in) {
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &proc_num_);

  SetTypeOfTask(GetStaticTypeOfTask());

  if (proc_rank_ == 0) {
    GetInput() = in;
  }

  GetOutput() = std::vector<std::pair<double, double>>();
}

bool PerepelkinIConvexHullGrahamScanALL::ValidationImpl() {
  return GetOutput().empty();
}

bool PerepelkinIConvexHullGrahamScanALL::PreProcessingImpl() {
  return true;
}

bool PerepelkinIConvexHullGrahamScanALL::RunImpl() {
  MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_POINT_);
  MPI_Type_commit(&MPI_POINT_);

  // [1] Broadcast original and padded sizes
  size_t original_size = 0;
  size_t padded_size = 0;
  BcastSizes(original_size, padded_size);

  // Handle edge cases
  if (original_size <= 1) {
    if (proc_rank_ == 0) {
      GetOutput() = GetInput();
    }
    BcastOutput();
    MPI_Type_free(&MPI_POINT_);
    return true;
  }

  tbb::global_control gc(tbb::global_control::max_allowed_parallelism, ppc::util::GetNumThreads());

  // [2] Distribute data
  std::vector<std::pair<double, double>> local_data;
  DistributeData(original_size, padded_size, local_data);

  // [3] Local find pivot

  // Remove fake points
  auto result = std::ranges::remove_if(local_data, [](const auto &p) { return p.first > 1e17; });
  local_data.erase(result.begin(), result.end());

  std::pair<double, double> local_pivot = {1e18, 1e18};
  FindPivotParallel(local_data, local_pivot);

  // [4] Find global pivot
  std::pair<double, double> global_pivot{};
  FindGlobalPivot(local_data, local_pivot, global_pivot);

  // [5] Local parallel sorting
  ParallelSort(local_data, global_pivot);

  // [6] Gather sorted blocks
  std::vector<int> real_counts(proc_num_);
  std::vector<int> real_displs(proc_num_);
  std::vector<std::pair<double, double>> gathered;
  GatherSortedBlocks(local_data, real_counts, real_displs, gathered);

  // [7] Merge sorted blocks (parallel)
  std::vector<std::pair<double, double>> sorted;
  if (proc_rank_ == 0) {
    sorted = MergeSortedBlocksParallel(gathered, real_counts, real_displs, global_pivot);
  }

  // [8] Sequential hull construction
  if (proc_rank_ == 0) {
    std::vector<std::pair<double, double>> hull;
    HullConstruction(hull, sorted, global_pivot);
    GetOutput() = std::move(hull);
  }

  // [9] Bcast output to all processes
  BcastOutput();

  MPI_Type_free(&MPI_POINT_);
  return true;
}

//
// Transfer data
//
void PerepelkinIConvexHullGrahamScanALL::BcastSizes(size_t &original_size, size_t &padded_size) {
  if (proc_rank_ == 0) {
    original_size = GetInput().size();
    const size_t remainder = original_size % proc_num_;
    padded_size = original_size + (remainder == 0 ? 0 : (proc_num_ - remainder));
  }

  MPI_Bcast(&original_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&padded_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
}

void PerepelkinIConvexHullGrahamScanALL::BcastOutput() {
  size_t size = 0;

  if (proc_rank_ == 0) {
    size = GetOutput().size();
  }

  MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

  GetOutput().resize(size);

  if (size > 0) {
    MPI_Bcast(GetOutput().data(), static_cast<int>(size), MPI_POINT_, 0, MPI_COMM_WORLD);
  }
}

void PerepelkinIConvexHullGrahamScanALL::DistributeData(const size_t &original_size, const size_t &padded_size,
                                                        std::vector<std::pair<double, double>> &local_data) {
  std::vector<std::pair<double, double>> padded_input;
  if (proc_rank_ == 0) {
    padded_input = GetInput();
    if (padded_size > original_size) {
      padded_input.resize(padded_size, std::make_pair(1e18, 1e18));
    }
  }

  std::vector<int> counts;
  std::vector<int> displs;
  const int base_size = static_cast<int>(padded_size / proc_num_);

  counts.resize(proc_num_);
  displs.resize(proc_num_);

  for (int i = 0, offset = 0; i < proc_num_; i++) {
    counts[i] = base_size;
    displs[i] = offset;
    offset += base_size;
  }

  const int local_size = counts[proc_rank_];
  local_data.resize(local_size);

  MPI_Scatterv(padded_input.data(), counts.data(), displs.data(), MPI_POINT_, local_data.data(), local_size, MPI_POINT_,
               0, MPI_COMM_WORLD);
}

void PerepelkinIConvexHullGrahamScanALL::FindGlobalPivot(std::vector<std::pair<double, double>> &local_data,
                                                         const std::pair<double, double> &local_pivot,
                                                         std::pair<double, double> &global_pivot) {
  std::vector<std::pair<double, double>> pivots;
  if (proc_rank_ == 0) {
    pivots.resize(proc_num_);
  }

  MPI_Gather(&local_pivot, 1, MPI_POINT_, pivots.data(), 1, MPI_POINT_, 0, MPI_COMM_WORLD);

  // [4.2] Find global pivot

  if (proc_rank_ == 0) {
    global_pivot = pivots[0];

    for (int i = 1; i < proc_num_; i++) {
      if (IsBetterPivot(pivots[i], global_pivot)) {
        global_pivot = pivots[i];
      }
    }
  }

  MPI_Bcast(&global_pivot, 1, MPI_POINT_, 0, MPI_COMM_WORLD);

  // Remove global pivot
  int owner_rank = -1;
  if (local_pivot == global_pivot) {
    owner_rank = proc_rank_;
  }
  int pivot_rank = -1;
  MPI_Allreduce(&owner_rank, &pivot_rank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (proc_rank_ == pivot_rank) {
    auto it = std::ranges::find(local_data, global_pivot);
    if (it != local_data.end()) {
      local_data.erase(it);
    }
  }
}

void PerepelkinIConvexHullGrahamScanALL::GatherSortedBlocks(std::vector<std::pair<double, double>> &local_data,
                                                            std::vector<int> &real_counts,
                                                            std::vector<int> &real_displs,
                                                            std::vector<std::pair<double, double>> &gathered) {
  int local_size = static_cast<int>(local_data.size());

  MPI_Gather(&local_size, 1, MPI_INT, real_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (proc_rank_ == 0) {
    int offset = 0;
    for (int i = 0; i < proc_num_; i++) {
      real_displs[i] = offset;
      offset += real_counts[i];
    }
    gathered.resize(offset);
  }

  MPI_Gatherv(local_data.data(), local_size, MPI_POINT_, gathered.data(), real_counts.data(), real_displs.data(),
              MPI_POINT_, 0, MPI_COMM_WORLD);
}

//
// Merge blocks
//
std::vector<std::pair<double, double>> PerepelkinIConvexHullGrahamScanALL::MergeSortedBlocksParallel(
    const std::vector<std::pair<double, double>> &gathered, const std::vector<int> &counts,
    const std::vector<int> &displs, const std::pair<double, double> &pivot) {
  std::vector<std::vector<std::pair<double, double>>> blocks;
  blocks.reserve(proc_num_);

  for (int i = 0; i < proc_num_; i++) {
    if (counts[i] > 0) {
      blocks.emplace_back(gathered.begin() + displs[i], gathered.begin() + displs[i] + counts[i]);
    } else {
      blocks.emplace_back();
    }
  }

  return MergeBlocksRange(blocks, 0, static_cast<int>(blocks.size()), pivot);
}

std::vector<std::pair<double, double>> PerepelkinIConvexHullGrahamScanALL::MergeBlocksRange(
    const std::vector<std::vector<std::pair<double, double>>> &blocks, int left, int right,
    const std::pair<double, double> &pivot) {
  if (right - left <= 0) {
    return {};
  }

  if (right - left == 1) {
    return blocks[left];
  }

  const int mid = left + ((right - left) / 2);

  std::vector<std::pair<double, double>> merged_left;
  std::vector<std::pair<double, double>> merged_right;

  tbb::parallel_invoke([&]() { merged_left = MergeBlocksRange(blocks, left, mid, pivot); },
                       [&]() { merged_right = MergeBlocksRange(blocks, mid, right, pivot); });

  return MergeTwoBlocks(merged_left, merged_right, pivot);
}

std::vector<std::pair<double, double>> PerepelkinIConvexHullGrahamScanALL::MergeTwoBlocks(
    const std::vector<std::pair<double, double>> &left, const std::vector<std::pair<double, double>> &right,
    const std::pair<double, double> &pivot) {
  std::vector<std::pair<double, double>> result;
  result.reserve(left.size() + right.size());

  size_t i = 0;
  size_t j = 0;

  while (i < left.size() && j < right.size()) {
    if (AngleCmp(left[i], right[j], pivot)) {
      result.push_back(left[i]);
      i++;
    } else {
      result.push_back(right[j]);
      j++;
    }
  }

  while (i < left.size()) {
    result.push_back(left[i]);
    i++;
  }

  while (j < right.size()) {
    result.push_back(right[j]);
    j++;
  }

  return result;
}

//
// Threads parallelization
//
void PerepelkinIConvexHullGrahamScanALL::FindPivotParallel(const std::vector<std::pair<double, double>> &pts,
                                                           std::pair<double, double> &local_pivot) {
  if (pts.empty()) {
    return;
  }

  auto better = [&](size_t a, size_t b) {
    if (pts[b].second < pts[a].second || (pts[b].second == pts[a].second && pts[b].first < pts[a].first)) {
      return b;
    }
    return a;
  };

  size_t pivot_idx = tbb::parallel_reduce(tbb::blocked_range<size_t>(1, pts.size()), static_cast<size_t>(0),
                                          [&](const tbb::blocked_range<size_t> &r, size_t local_idx) {
    for (size_t i = r.begin(); i != r.end(); i++) {
      local_idx = better(local_idx, i);
    }
    return local_idx;
  }, [&](size_t a, size_t b) { return better(a, b); });

  local_pivot = pts[pivot_idx];
}

void PerepelkinIConvexHullGrahamScanALL::ParallelSort(std::vector<std::pair<double, double>> &data,
                                                      const std::pair<double, double> &pivot) {
  tbb::parallel_sort(data.begin(), data.end(), [&](const auto &a, const auto &b) { return AngleCmp(a, b, pivot); });
}

//
// Sequential
//
void PerepelkinIConvexHullGrahamScanALL::HullConstruction(std::vector<std::pair<double, double>> &hull,
                                                          const std::vector<std::pair<double, double>> &pts,
                                                          const std::pair<double, double> &pivot) {
  hull.clear();
  hull.push_back(pivot);

  if (pts.empty()) {
    return;
  }

  hull.push_back(pts[0]);

  for (size_t i = 1; i < pts.size(); i++) {
    while (hull.size() >= 2 && Orientation(hull[hull.size() - 2], hull[hull.size() - 1], pts[i]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(pts[i]);
  }
}

//
// Helpers
//
bool PerepelkinIConvexHullGrahamScanALL::IsBetterPivot(const std::pair<double, double> &a,
                                                       const std::pair<double, double> &b) {
  return (a.second < b.second) || (a.second == b.second && a.first < b.first);
}

double PerepelkinIConvexHullGrahamScanALL::Orientation(const std::pair<double, double> &p,
                                                       const std::pair<double, double> &q,
                                                       const std::pair<double, double> &r) {
  double val = ((q.first - p.first) * (r.second - p.second)) - ((q.second - p.second) * (r.first - p.first));

  if (std::abs(val) < 1e-9) {
    return 0.0;
  }

  return val;
}

bool PerepelkinIConvexHullGrahamScanALL::AngleCmp(const std::pair<double, double> &a,
                                                  const std::pair<double, double> &b,
                                                  const std::pair<double, double> &pivot) {
  double dx1 = a.first - pivot.first;
  double dy1 = a.second - pivot.second;
  double dx2 = b.first - pivot.first;
  double dy2 = b.second - pivot.second;

  double cross = (dx1 * dy2) - (dy1 * dx2);

  if (std::abs(cross) < 1e-9) {
    double dist1 = (dx1 * dx1) + (dy1 * dy1);
    double dist2 = (dx2 * dx2) + (dy2 * dy2);
    return dist1 < dist2;
  }

  return cross > 0;
}

bool PerepelkinIConvexHullGrahamScanALL::PostProcessingImpl() {
  return true;
}

}  // namespace perepelkin_i_convex_hull_graham_scan
