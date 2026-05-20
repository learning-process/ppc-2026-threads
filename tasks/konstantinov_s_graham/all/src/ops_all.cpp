#include "konstantinov_s_graham/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <utility>
#include <vector>

#include "konstantinov_s_graham/common/include/common.hpp"
#include "util/include/util.hpp"

namespace konstantinov_s_graham {

KonstantinovAGrahamALL::KonstantinovAGrahamALL(const InType &in) {
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &proc_num_);

  SetTypeOfTask(GetStaticTypeOfTask());

  if (proc_rank_ == 0) {
    GetInput() = in;
  }

  GetOutput() = OutType();
}

bool KonstantinovAGrahamALL::ValidationImpl() {
  return GetInput().first.size() == GetInput().second.size();
}

bool KonstantinovAGrahamALL::PreProcessingImpl() {
  return true;
}

void KonstantinovAGrahamALL::RemoveDuplicates(std::vector<double> &xs, std::vector<double> &ys) {
  std::vector<std::pair<double, double>> pts;
  pts.reserve(xs.size());

  for (size_t i = 0; i < xs.size(); ++i) {
    pts.emplace_back(xs[i], ys[i]);
  }

  std::ranges::sort(pts, [](const auto &lhs, const auto &rhs) {
    if (std::abs(lhs.first - rhs.first) > kKEps) {
      return lhs.first < rhs.first;
    }
    return lhs.second < rhs.second;
  });

  const auto unique_end = std::ranges::unique(pts, [](const auto &lhs, const auto &rhs) {
    return std::abs(lhs.first - rhs.first) < kKEps && std::abs(lhs.second - rhs.second) < kKEps;
  });

  pts.erase(unique_end.begin(), pts.end());

  xs.resize(pts.size());
  ys.resize(pts.size());

  for (size_t i = 0; i < pts.size(); ++i) {
    xs[i] = pts[i].first;
    ys[i] = pts[i].second;
  }
}

bool KonstantinovAGrahamALL::IsLowerAnchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t lhs,
                                           size_t rhs) {
  if (ys[lhs] < ys[rhs] - kKEps) {
    return true;
  }

  if (std::abs(ys[lhs] - ys[rhs]) < kKEps && xs[lhs] < xs[rhs] - kKEps) {
    return true;
  }

  return false;
}

size_t KonstantinovAGrahamALL::FindAnchorIndex(const std::vector<double> &xs, const std::vector<double> &ys) {
  if (xs.empty()) {
    return 0;
  }

  return tbb::parallel_reduce(tbb::blocked_range<size_t>(1, xs.size()), size_t{0},
                              [&xs, &ys](const tbb::blocked_range<size_t> &range, size_t local_idx) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      if (IsLowerAnchor(xs, ys, i, local_idx)) {
        local_idx = i;
      }
    }
    return local_idx;
  }, [&xs, &ys](size_t left, size_t right) { return IsLowerAnchor(xs, ys, left, right) ? left : right; });
}

double KonstantinovAGrahamALL::Dist2(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j) {
  const double dx = xs[j] - xs[i];
  const double dy = ys[j] - ys[i];
  return (dx * dx) + (dy * dy);
}

double KonstantinovAGrahamALL::CrossVal(const std::vector<double> &xs, const std::vector<double> &ys, size_t i,
                                        size_t j, size_t k) {
  const double ax = xs[j] - xs[i];
  const double ay = ys[j] - ys[i];
  const double bx = xs[k] - xs[i];
  const double by = ys[k] - ys[i];
  return (ax * by) - (ay * bx);
}

void KonstantinovAGrahamALL::FillIndexRange(std::vector<size_t> &idxs, size_t begin, size_t end, size_t anchor_idx) {
  for (size_t i = begin; i < end; ++i) {
    if (i < anchor_idx) {
      idxs[i] = i;
    } else if (i > anchor_idx) {
      idxs[i - 1] = i;
    }
  }
}

void KonstantinovAGrahamALL::FillIndicesParallel(std::vector<size_t> &idxs, size_t point_count, size_t anchor_idx) {
  const auto thread_count = static_cast<size_t>(ppc::util::GetNumThreads());

  if (thread_count <= 1 || point_count < (thread_count * 2)) {
    FillIndexRange(idxs, 0, point_count, anchor_idx);
    return;
  }

  std::vector<tbb::blocked_range<size_t>> ranges;
  ranges.reserve(thread_count);

  const size_t block_size = (point_count + thread_count - 1) / thread_count;
  for (size_t thread_idx = 0; thread_idx < thread_count; ++thread_idx) {
    const size_t begin = thread_idx * block_size;
    const size_t end = std::min(point_count, begin + block_size);
    if (begin < end) {
      ranges.emplace_back(begin, end);
    }
  }

  tbb::parallel_for(tbb::blocked_range<size_t>(0, ranges.size()), [&](const tbb::blocked_range<size_t> &outer_range) {
    for (size_t i = outer_range.begin(); i < outer_range.end(); ++i) {
      FillIndexRange(idxs, ranges[i].begin(), ranges[i].end(), anchor_idx);
    }
  });
}

std::vector<size_t> KonstantinovAGrahamALL::CollectAndSortIndices(const std::vector<double> &xs,
                                                                  const std::vector<double> &ys, size_t anchor_idx) {
  std::vector<size_t> idxs(xs.size() - 1);

  FillIndicesParallel(idxs, xs.size(), anchor_idx);

  tbb::parallel_sort(idxs.begin(), idxs.end(), [&xs, &ys, anchor_idx](size_t lhs, size_t rhs) {
    const double cross = CrossVal(xs, ys, anchor_idx, lhs, rhs);

    if (std::abs(cross) < kKEps) {
      return Dist2(xs, ys, anchor_idx, lhs) < Dist2(xs, ys, anchor_idx, rhs);
    }

    return cross > 0;
  });

  return idxs;
}

bool KonstantinovAGrahamALL::CheckCollinearRange(const std::vector<double> &xs, const std::vector<double> &ys,
                                                 size_t anchor_idx, const std::vector<size_t> &sorted_idxs,
                                                 size_t begin, size_t end) {
  for (size_t i = begin; i < end; ++i) {
    if (std::abs(CrossVal(xs, ys, anchor_idx, sorted_idxs[0], sorted_idxs[i])) > kKEps) {
      return false;
    }
  }

  return true;
}

bool KonstantinovAGrahamALL::AllCollinearWithAnchor(const std::vector<double> &xs, const std::vector<double> &ys,
                                                    size_t anchor_idx, const std::vector<size_t> &sorted_idxs) {
  if (sorted_idxs.size() < 2) {
    return true;
  }

  const auto thread_count = static_cast<size_t>(ppc::util::GetNumThreads());

  if (thread_count <= 1 || sorted_idxs.size() < (thread_count * 2)) {
    return CheckCollinearRange(xs, ys, anchor_idx, sorted_idxs, 1, sorted_idxs.size());
  }

  std::vector<tbb::blocked_range<size_t>> ranges;
  ranges.reserve(thread_count);

  const size_t block_size = (sorted_idxs.size() + thread_count - 1) / thread_count;
  for (size_t thread_idx = 0; thread_idx < thread_count; ++thread_idx) {
    const size_t begin = std::max<size_t>(1, thread_idx * block_size);
    const size_t end = std::min(sorted_idxs.size(), begin + block_size);
    if (begin < end) {
      ranges.emplace_back(begin, end);
    }
  }

  return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, ranges.size()), true,
                              [&](const tbb::blocked_range<size_t> &outer_range, bool local_ok) {
    for (size_t i = outer_range.begin(); i < outer_range.end() && local_ok; ++i) {
      if (!CheckCollinearRange(xs, ys, anchor_idx, sorted_idxs, ranges[i].begin(), ranges[i].end())) {
        local_ok = false;
      }
    }
    return local_ok;
  }, [](bool lhs, bool rhs) { return lhs && rhs; });
}

std::vector<std::pair<double, double>> KonstantinovAGrahamALL::BuildHullFromSorted(
    const std::vector<double> &xs, const std::vector<double> &ys, size_t anchor_idx,
    const std::vector<size_t> &sorted_idxs) {
  std::vector<size_t> stack;
  stack.reserve(sorted_idxs.size() + 1);
  stack.push_back(anchor_idx);

  if (!sorted_idxs.empty()) {
    stack.push_back(sorted_idxs[0]);
  }

  for (size_t i = 1; i < sorted_idxs.size(); ++i) {
    const size_t cur = sorted_idxs[i];

    while (stack.size() >= 2) {
      const size_t q = stack.back();
      const size_t p = stack[stack.size() - 2];
      const double cross = CrossVal(xs, ys, p, q, cur);

      if (cross <= kKEps) {
        stack.pop_back();
      } else {
        break;
      }
    }

    stack.push_back(cur);
  }

  std::vector<std::pair<double, double>> hull;
  hull.reserve(stack.size());

  for (size_t id : stack) {
    hull.emplace_back(xs[id], ys[id]);
  }

  return hull;
}

std::vector<std::pair<double, double>> KonstantinovAGrahamALL::BuildHullFromCoords(const std::vector<double> &xs,
                                                                                   const std::vector<double> &ys) {
  std::vector<double> local_xs = xs;
  std::vector<double> local_ys = ys;

  RemoveDuplicates(local_xs, local_ys);

  if (local_xs.size() != local_ys.size() || local_xs.empty()) {
    return {};
  }

  if (local_xs.size() < 3) {
    std::vector<std::pair<double, double>> out;
    out.reserve(local_xs.size());

    for (size_t i = 0; i < local_xs.size(); ++i) {
      out.emplace_back(local_xs[i], local_ys[i]);
    }

    return out;
  }

  const size_t anchor_idx = FindAnchorIndex(local_xs, local_ys);
  const std::vector<size_t> sorted_idxs = CollectAndSortIndices(local_xs, local_ys, anchor_idx);

  if (sorted_idxs.empty()) {
    return {{local_xs[anchor_idx], local_ys[anchor_idx]}};
  }

  if (AllCollinearWithAnchor(local_xs, local_ys, anchor_idx, sorted_idxs)) {
    const size_t far_idx = sorted_idxs.back();
    return {{local_xs[anchor_idx], local_ys[anchor_idx]}, {local_xs[far_idx], local_ys[far_idx]}};
  }

  return BuildHullFromSorted(local_xs, local_ys, anchor_idx, sorted_idxs);
}

void KonstantinovAGrahamALL::ScatterInput(size_t total_size, std::vector<double> &local_xs,
                                          std::vector<double> &local_ys) {
  std::vector<int> counts(proc_num_, 0);
  std::vector<int> displs(proc_num_, 0);

  const size_t base = total_size / static_cast<size_t>(proc_num_);
  const size_t remainder = total_size % static_cast<size_t>(proc_num_);

  size_t offset = 0;
  for (int rank = 0; rank < proc_num_; ++rank) {
    const size_t amount = base + (std::cmp_less(rank, remainder) ? 1U : 0U);
    counts[rank] = static_cast<int>(amount);
    displs[rank] = static_cast<int>(offset);
    offset += amount;
  }

  const int local_size = counts[proc_rank_];
  local_xs.resize(static_cast<size_t>(local_size));
  local_ys.resize(static_cast<size_t>(local_size));

  const double *send_xs = proc_rank_ == 0 ? GetInput().first.data() : nullptr;
  const double *send_ys = proc_rank_ == 0 ? GetInput().second.data() : nullptr;

  MPI_Scatterv(send_xs, counts.data(), displs.data(), MPI_DOUBLE, local_xs.data(), local_size, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
  MPI_Scatterv(send_ys, counts.data(), displs.data(), MPI_DOUBLE, local_ys.data(), local_size, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
}

void KonstantinovAGrahamALL::GatherLocalHull(const std::vector<std::pair<double, double>> &local_hull,
                                             std::vector<double> &gathered_xs, std::vector<double> &gathered_ys) const {
  const auto local_size = static_cast<int>(local_hull.size());

  std::vector<int> counts(proc_num_, 0);
  MPI_Gather(&local_size, 1, MPI_INT, proc_rank_ == 0 ? counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displs(proc_num_, 0);
  int total_size = 0;

  if (proc_rank_ == 0) {
    for (int rank = 0; rank < proc_num_; ++rank) {
      displs[rank] = total_size;
      total_size += counts[rank];
    }

    gathered_xs.resize(static_cast<size_t>(total_size));
    gathered_ys.resize(static_cast<size_t>(total_size));
  }

  std::vector<double> local_xs(static_cast<size_t>(local_size));
  std::vector<double> local_ys(static_cast<size_t>(local_size));

  for (int i = 0; i < local_size; ++i) {
    local_xs[static_cast<size_t>(i)] = local_hull[static_cast<size_t>(i)].first;
    local_ys[static_cast<size_t>(i)] = local_hull[static_cast<size_t>(i)].second;
  }

  MPI_Gatherv(local_xs.data(), local_size, MPI_DOUBLE, proc_rank_ == 0 ? gathered_xs.data() : nullptr,
              proc_rank_ == 0 ? counts.data() : nullptr, proc_rank_ == 0 ? displs.data() : nullptr, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  MPI_Gatherv(local_ys.data(), local_size, MPI_DOUBLE, proc_rank_ == 0 ? gathered_ys.data() : nullptr,
              proc_rank_ == 0 ? counts.data() : nullptr, proc_rank_ == 0 ? displs.data() : nullptr, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
}

void KonstantinovAGrahamALL::BroadcastOutput() {
  std::uint64_t size_u64 = 0;

  if (proc_rank_ == 0) {
    size_u64 = static_cast<std::uint64_t>(GetOutput().size());
  }

  MPI_Bcast(&size_u64, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  const auto size = static_cast<size_t>(size_u64);

  std::vector<double> xs(size);
  std::vector<double> ys(size);

  if (proc_rank_ == 0) {
    for (size_t i = 0; i < size; ++i) {
      xs[i] = GetOutput()[i].first;
      ys[i] = GetOutput()[i].second;
    }
  }

  if (size > 0) {
    MPI_Bcast(xs.data(), static_cast<int>(size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ys.data(), static_cast<int>(size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  if (proc_rank_ != 0) {
    GetOutput().resize(size);
    for (size_t i = 0; i < size; ++i) {
      GetOutput()[i] = {xs[i], ys[i]};
    }
  }
}

bool KonstantinovAGrahamALL::RunImpl() {
  const tbb::global_control gc(tbb::global_control::max_allowed_parallelism, ppc::util::GetNumThreads());

  std::uint64_t total_size_u64 = 0;

  if (proc_rank_ == 0) {
    total_size_u64 = static_cast<std::uint64_t>(GetInput().first.size());
  }

  MPI_Bcast(&total_size_u64, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  const auto total_size = static_cast<size_t>(total_size_u64);

  if (total_size == 0) {
    if (proc_rank_ == 0) {
      GetOutput() = {};
    }
    BroadcastOutput();
    return true;
  }

  if (total_size < 3 || proc_num_ == 1) {
    if (proc_rank_ == 0) {
      GetOutput() = BuildHullFromCoords(GetInput().first, GetInput().second);
    }
    BroadcastOutput();
    return true;
  }

  std::vector<double> local_xs;
  std::vector<double> local_ys;
  ScatterInput(total_size, local_xs, local_ys);

  const std::vector<std::pair<double, double>> local_hull = BuildHullFromCoords(local_xs, local_ys);

  std::vector<double> gathered_xs;
  std::vector<double> gathered_ys;
  GatherLocalHull(local_hull, gathered_xs, gathered_ys);

  if (proc_rank_ == 0) {
    GetOutput() = BuildHullFromCoords(gathered_xs, gathered_ys);
  }

  BroadcastOutput();
  return true;
}

bool KonstantinovAGrahamALL::PostProcessingImpl() {
  return true;
}

}  // namespace konstantinov_s_graham
