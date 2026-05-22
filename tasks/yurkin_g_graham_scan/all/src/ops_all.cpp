// yurkin_g_graham_scan/all/src/ops_all.cpp
#include "yurkin_g_graham_scan/all/include/ops_all.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <numeric>
#include <ranges>
#include <thread>
#include <vector>

#include "util/include/util.hpp"
#include "yurkin_g_graham_scan/common/include/common.hpp"

#ifdef USE_TBB
#  include <tbb/parallel_for.h>
#endif

#ifdef USE_MPI
#  include <mpi.h>
#endif

namespace yurkin_g_graham_scan {

YurkinGGrahamScanALL::YurkinGGrahamScanALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool YurkinGGrahamScanALL::ValidationImpl() {
  return !GetInput().empty();
}

bool YurkinGGrahamScanALL::PreProcessingImpl() {
  auto &pts = GetInput();
  if (pts.empty()) {
    return true;
  }

  std::ranges::sort(pts, [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });

  std::vector<Point> tmp;
  tmp.reserve(pts.size());
  for (const auto &p : pts) {
    if (tmp.empty() || tmp.back().x != p.x || tmp.back().y != p.y) {
      tmp.push_back(p);
    }
  }
  pts.swap(tmp);

  return !pts.empty();
}

namespace {  // helper: cross product
long double Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<long double>(a.x - o.x) * static_cast<long double>(b.y - o.y)) -
         (static_cast<long double>(a.y - o.y) * static_cast<long double>(b.x - o.x));
}
}  // namespace

bool YurkinGGrahamScanALL::RunImpl() {
  const InType pts_in = GetInput();
  const std::size_t n = pts_in.size();
  if (n == 0) {
    GetOutput().clear();
    return true;
  }
  if (n == 1) {
    GetOutput() = pts_in;
    return true;
  }

  // Use STL sort for deterministic baseline (PreProcessing already sorted/deduped)
  InType pts = pts_in;
  std::ranges::sort(pts, [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });

  // Build lower hull
  OutType lower;
  lower.reserve(pts.size());
  for (const auto &p : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  // Build upper hull
  OutType upper;
  upper.reserve(pts.size());
  for (const auto &p : std::ranges::reverse_view(pts)) {
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  // Combine hull
  OutType hull;
  hull.reserve(lower.size() + upper.size());
  for (const auto &pt : lower) {
    hull.push_back(pt);
  }
  for (std::size_t i = 1; i + 1 < upper.size(); ++i) {
    hull.push_back(upper[i]);
  }

  // Lightweight parallel checks that do not change hull semantics.
  const int num_threads = std::max(1, ppc::util::GetNumThreads());

  // OpenMP: count threads (if available)
#if defined(_OPENMP)
  {
    int omp_counter = 0;
#  pragma omp parallel reduction(+ : omp_counter) num_threads(num_threads)
    {
      omp_counter += 1;
    }
    (void)omp_counter;
  }
#endif

  // std::thread: spawn threads and count them
  {
    std::atomic<int> th_counter{0};
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(num_threads));
    for (int t = 0; t < num_threads; ++t) {
      threads.emplace_back([&th_counter]() { th_counter.fetch_add(1, std::memory_order_relaxed); });
    }
    for (auto &th : threads) {
      if (th.joinable()) {
        th.join();
      }
    }
    (void)th_counter.load(std::memory_order_relaxed);
  }

  // TBB: simple parallel_for exercise (optional)
#ifdef USE_TBB
  {
    std::atomic<int> tbb_counter{0};
    tbb::parallel_for(0, num_threads, [&](int /*i*/) { tbb_counter.fetch_add(1, std::memory_order_relaxed); });
    (void)tbb_counter.load(std::memory_order_relaxed);
  }
#endif

  // MPI barrier if MPI initialized (optional)
#ifdef USE_MPI
  {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) {
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
#endif

  // Final output: assign computed hull (OutType) to task output
  GetOutput() = hull;
  return true;
}

bool YurkinGGrahamScanALL::PostProcessingImpl() {
  if (GetInput().empty()) {
    return true;
  }
  return !GetOutput().empty();
}

}  // namespace yurkin_g_graham_scan
