#include "dergachev_a_graham_scan/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "dergachev_a_graham_scan/common/include/common.hpp"
#include "util/include/util.hpp"

namespace dergachev_a_graham_scan {

namespace {

using Pt = std::pair<double, double>;

constexpr double kPaddingMarker = 1e18;

const double kPi = std::acos(-1.0);

struct Slice {
  const Pt *begin;
  const Pt *end;
};

struct WorkBuffers {
  std::vector<Pt> padded_input;
  std::vector<Pt> local_data;
  std::vector<Pt> gathered_pivots;
  std::vector<int> counts;
  std::vector<int> displs;
  std::vector<int> recv_counts;
  std::vector<int> recv_displs;
  std::vector<Pt> gathered;
  std::vector<Pt> sorted;
  std::vector<Pt> merge_temp;
};

MPI_Datatype GetMpiPointType() {
  static MPI_Datatype mpi_point = MPI_DATATYPE_NULL;
  static bool initialized = false;
  if (!initialized) {
    MPI_Type_contiguous(2, MPI_DOUBLE, &mpi_point);
    MPI_Type_commit(&mpi_point);
    initialized = true;
  }
  return mpi_point;
}

bool IsLowerLeft(const Pt &a, const Pt &b) {
  return a.second < b.second || (a.second == b.second && a.first < b.first);
}

double CrossProduct(const Pt &o, const Pt &a, const Pt &b) {
  return ((a.first - o.first) * (b.second - o.second)) - ((a.second - o.second) * (b.first - o.first));
}

double DistSquared(const Pt &a, const Pt &b) {
  double dx = a.first - b.first;
  double dy = a.second - b.second;
  return (dx * dx) + (dy * dy);
}

bool AngleLessSeq(const Pt &a, const Pt &b, const Pt &pivot) {
  const double cross = CrossProduct(pivot, a, b);
  return cross > 0.0 || (cross == 0.0 && DistSquared(pivot, a) < DistSquared(pivot, b));
}

int FindLocalPivotIndex(const std::vector<Pt> &pts) {
  if (pts.empty()) {
    return 0;
  }

  const int num_threads = ppc::util::GetNumThreads();
  const int n = static_cast<int>(pts.size());
  if (n < num_threads * 2) {
    int best = 0;
    for (int i = 1; i < n; i++) {
      if (IsLowerLeft(pts[i], pts[best])) {
        best = i;
      }
    }
    return best;
  }

  std::vector<int> local_best(num_threads);
#pragma omp parallel num_threads(num_threads) default(none) shared(pts, n, local_best, num_threads)
  {
    const int tid = omp_get_thread_num();
    const int lo = (tid * n) / num_threads;
    const int hi = ((tid + 1) * n) / num_threads;
    int best = lo;
    for (int i = lo + 1; i < hi; i++) {
      if (IsLowerLeft(pts[i], pts[best])) {
        best = i;
      }
    }
    local_best[tid] = best;
  }

  int best = local_best[0];
  for (int thread_idx = 1; thread_idx < num_threads; thread_idx++) {
    if (IsLowerLeft(pts[local_best[thread_idx]], pts[best])) {
      best = local_best[thread_idx];
    }
  }
  return best;
}

void ParallelSortRange(std::vector<Pt>::iterator begin, std::vector<Pt>::iterator end, const Pt &pivot) {
  const int n = static_cast<int>(end - begin);
  const int num_threads = ppc::util::GetNumThreads();
  if (n <= 1 || num_threads <= 1) {
    std::sort(begin, end, [&](const Pt &a, const Pt &b) { return AngleLessSeq(a, b, pivot); });
    return;
  }

  auto cmp = [&](const Pt &a, const Pt &b) { return AngleLessSeq(a, b, pivot); };
  const int chunk = n / num_threads;
#pragma omp parallel num_threads(num_threads) default(none) shared(begin, n, chunk, cmp, num_threads)
  {
    const int tid = omp_get_thread_num();
    const int lo = tid * chunk;
    const int hi = (tid == num_threads - 1) ? n : (tid + 1) * chunk;
    std::sort(begin + lo, begin + hi, cmp);
  }

  int boundary = chunk;
  for (int tid = 1; tid < num_threads; tid++) {
    const int next = (tid == num_threads - 1) ? n : (tid + 1) * chunk;
    std::inplace_merge(begin, begin + boundary, begin + next, cmp);
    boundary = next;
  }
}

void BuildHullFromSorted(const std::vector<Pt> &sorted, const Pt &pivot, std::vector<Pt> &hull) {
  hull.clear();
  hull.push_back(pivot);
  for (const auto &p : sorted) {
    while (hull.size() >= 2 && CrossProduct(hull[hull.size() - 2], hull.back(), p) <= 0.0) {
      hull.pop_back();
    }
    hull.push_back(p);
  }
}

void MergeTwoSlices(Slice left, Slice right, const Pt &pivot, std::vector<Pt> &out) {
  out.clear();
  out.reserve(static_cast<size_t>((left.end - left.begin) + (right.end - right.begin)));

  const Pt *i = left.begin;
  const Pt *j = right.begin;
  while (i < left.end && j < right.end) {
    if (AngleLessSeq(*i, *j, pivot)) {
      out.push_back(*i++);
    } else {
      out.push_back(*j++);
    }
  }
  while (i < left.end) {
    out.push_back(*i++);
  }
  while (j < right.end) {
    out.push_back(*j++);
  }
}

Slice MakeSlice(const std::vector<Pt> &data, int displ, int count) {
  return {.begin = data.data() + displ, .end = data.data() + displ + count};
}

std::vector<Pt> MergeAllBlocks(const std::vector<Pt> &gathered, const std::vector<int> &displs,
                               const std::vector<int> &counts, int world_size, const Pt &pivot,
                               std::vector<Pt> &merge_temp) {
  const int first_displ = displs[0];
  const int first_count = counts[0];
  std::vector<Pt> merged(gathered.begin() + first_displ, gathered.begin() + first_displ + first_count);

  for (int block_idx = 1; block_idx < world_size; block_idx++) {
    const int displ = displs[static_cast<size_t>(block_idx)];
    const int count = counts[static_cast<size_t>(block_idx)];
    MergeTwoSlices(MakeSlice(merged, 0, static_cast<int>(merged.size())), MakeSlice(gathered, displ, count), pivot,
                   merge_temp);
    merged.swap(merge_temp);
  }

  return merged;
}

void RemovePaddingPoints(std::vector<Pt> &pts) {
  auto end_it = std::ranges::remove_if(pts, [](const Pt &p) { return p.first > kPaddingMarker / 10.0; });
  pts.erase(end_it.begin(), end_it.end());
}

void BcastHullSize(std::vector<Pt> &hull, int rank) {
  int hull_size = (rank == 0) ? static_cast<int>(hull.size()) : 0;
  MPI_Bcast(&hull_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  hull.resize(static_cast<size_t>(hull_size));
}

void ComputePaddedSizes(int world_size, int points_size, int &original_size, int &padded_size) {
  original_size = points_size;
  const int remainder = original_size % world_size;
  padded_size = original_size + ((remainder == 0) ? 0 : (world_size - remainder));
}

bool AllPointsSame(const std::vector<Pt> &points) {
  if (points.size() <= 1) {
    return false;
  }
  return std::all_of(points.begin() + 1, points.end(),
                     [&](const Pt &p) { return p.first == points[0].first && p.second == points[0].second; });
}

bool HandleTrivialCases(int rank, int original_size, int all_same, std::vector<Pt> &points, std::vector<Pt> &hull) {
  if (original_size <= 1) {
    if (rank == 0 && !points.empty()) {
      hull.push_back(points[0]);
    }
    BcastHullSize(hull, rank);
    return true;
  }

  if (all_same != 0) {
    if (rank == 0) {
      hull.push_back(points[0]);
    }
    BcastHullSize(hull, rank);
    return true;
  }

  return false;
}

void PrepareScatterLayout(WorkBuffers &bufs, int world_size, int block_size) {
  bufs.counts.assign(static_cast<size_t>(world_size), block_size);
  bufs.displs.resize(static_cast<size_t>(world_size));
  for (int i = 0, offset = 0; i < world_size; i++) {
    bufs.displs[static_cast<size_t>(i)] = offset;
    offset += block_size;
  }
}

void ScatterLocalPoints(WorkBuffers &bufs, int rank, int block_size, int padded_size, int original_size,
                        MPI_Datatype mpi_point, const std::vector<Pt> &points) {
  if (rank == 0) {
    bufs.padded_input = points;
    if (padded_size > original_size) {
      bufs.padded_input.resize(static_cast<size_t>(padded_size), Pt{kPaddingMarker, kPaddingMarker});
    }
  }

  bufs.local_data.resize(static_cast<size_t>(block_size));
  MPI_Scatterv(rank == 0 ? bufs.padded_input.data() : nullptr, bufs.counts.data(), bufs.displs.data(), mpi_point,
               bufs.local_data.data(), block_size, mpi_point, 0, MPI_COMM_WORLD);
  RemovePaddingPoints(bufs.local_data);
}

Pt FindGlobalPivot(WorkBuffers &bufs, int rank, int world_size, MPI_Datatype mpi_point, const Pt &local_pivot) {
  if (rank == 0) {
    bufs.gathered_pivots.resize(static_cast<size_t>(world_size));
  }
  MPI_Gather(&local_pivot, 1, mpi_point, rank == 0 ? bufs.gathered_pivots.data() : nullptr, 1, mpi_point, 0,
             MPI_COMM_WORLD);

  Pt global_pivot{kPaddingMarker, kPaddingMarker};
  if (rank == 0) {
    global_pivot = bufs.gathered_pivots[0];
    for (int i = 1; i < world_size; i++) {
      if (IsLowerLeft(bufs.gathered_pivots[static_cast<size_t>(i)], global_pivot)) {
        global_pivot = bufs.gathered_pivots[static_cast<size_t>(i)];
      }
    }
  }
  MPI_Bcast(&global_pivot, 1, mpi_point, 0, MPI_COMM_WORLD);
  return global_pivot;
}

void RemoveGlobalPivot(WorkBuffers &bufs, int rank, int pivot_owner, const Pt &global_pivot) {
  if (rank != pivot_owner) {
    return;
  }
  auto it = std::ranges::find(bufs.local_data, global_pivot);
  if (it != bufs.local_data.end()) {
    bufs.local_data.erase(it);
  }
}

void GatherSortedAndBuildHull(WorkBuffers &bufs, int rank, int world_size, MPI_Datatype mpi_point,
                              const Pt &global_pivot, std::vector<Pt> &hull) {
  const int local_size = static_cast<int>(bufs.local_data.size());
  MPI_Gather(&local_size, 1, MPI_INT, rank == 0 ? bufs.recv_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    int offset = 0;
    for (int i = 0; i < world_size; i++) {
      bufs.recv_displs[static_cast<size_t>(i)] = offset;
      offset += bufs.recv_counts[static_cast<size_t>(i)];
    }
    bufs.gathered.resize(static_cast<size_t>(offset));
  }

  MPI_Gatherv(bufs.local_data.data(), local_size, mpi_point, rank == 0 ? bufs.gathered.data() : nullptr,
              rank == 0 ? bufs.recv_counts.data() : nullptr, rank == 0 ? bufs.recv_displs.data() : nullptr, mpi_point,
              0, MPI_COMM_WORLD);

  if (rank != 0) {
    return;
  }

  bufs.sorted =
      MergeAllBlocks(bufs.gathered, bufs.recv_displs, bufs.recv_counts, world_size, global_pivot, bufs.merge_temp);
  BuildHullFromSorted(bufs.sorted, global_pivot, hull);
}

void RunDistributedHull(WorkBuffers &bufs, int rank, int world_size, int block_size, int padded_size, int original_size,
                        MPI_Datatype mpi_point, const std::vector<Pt> &points, std::vector<Pt> &hull) {
  PrepareScatterLayout(bufs, world_size, block_size);
  ScatterLocalPoints(bufs, rank, block_size, padded_size, original_size, mpi_point, points);

  const int local_pivot_idx = FindLocalPivotIndex(bufs.local_data);
  const Pt local_pivot = bufs.local_data.empty() ? Pt{kPaddingMarker, kPaddingMarker}
                                                 : bufs.local_data[static_cast<size_t>(local_pivot_idx)];

  const Pt global_pivot = FindGlobalPivot(bufs, rank, world_size, mpi_point, local_pivot);

  int owner_rank = (local_pivot == global_pivot) ? rank : -1;
  int pivot_owner = -1;
  MPI_Allreduce(&owner_rank, &pivot_owner, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  RemoveGlobalPivot(bufs, rank, pivot_owner, global_pivot);

  ParallelSortRange(bufs.local_data.begin(), bufs.local_data.end(), global_pivot);
  GatherSortedAndBuildHull(bufs, rank, world_size, mpi_point, global_pivot, hull);
  BcastHullSize(hull, rank);
}

}  // namespace

DergachevAGrahamScanALL::DergachevAGrahamScanALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool DergachevAGrahamScanALL::ValidationImpl() {
  return GetInput() >= 0;
}

bool DergachevAGrahamScanALL::PreProcessingImpl() {
  hull_.clear();
  int n_input = GetInput();
  if (n_input <= 0) {
    points_.clear();
    return true;
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    points_.clear();
    return true;
  }

  points_.resize(n_input);
  double step = (2.0 * kPi) / n_input;
  auto *pts_data = points_.data();
#pragma omp parallel for default(none) shared(pts_data, step, n_input) num_threads(ppc::util::GetNumThreads())
  for (int i = 0; i < n_input; i++) {
    pts_data[i] = {std::cos(step * i), std::sin(step * i)};
  }
  if (n_input > 3) {
    points_.emplace_back(0.0, 0.0);
  }
  return true;
}

bool DergachevAGrahamScanALL::RunImpl() {
  hull_.clear();

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Datatype mpi_point = GetMpiPointType();

  int original_size = 0;
  int padded_size = 0;
  if (rank == 0) {
    ComputePaddedSizes(world_size, static_cast<int>(points_.size()), original_size, padded_size);
  }
  MPI_Bcast(&original_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&padded_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int all_same = 0;
  if (rank == 0 && original_size > 1) {
    all_same = AllPointsSame(points_) ? 1 : 0;
  }
  MPI_Bcast(&all_same, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (HandleTrivialCases(rank, original_size, all_same, points_, hull_)) {
    return true;
  }

  thread_local WorkBuffers bufs;
  bufs.recv_counts.resize(static_cast<size_t>(world_size));
  bufs.recv_displs.resize(static_cast<size_t>(world_size));

  const int block_size = padded_size / world_size;
  RunDistributedHull(bufs, rank, world_size, block_size, padded_size, original_size, mpi_point, points_, hull_);
  return true;
}

bool DergachevAGrahamScanALL::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan
