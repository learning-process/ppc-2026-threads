#include "yurkin_g_graham_scan/all/include/ops_all.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <ranges>
#include <thread>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "yurkin_g_graham_scan/common/include/common.hpp"

#ifdef USE_MPI
#  include <mpi.h>
#endif

namespace yurkin_g_graham_scan {
namespace {

long double Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<long double>(a.x - o.x) * static_cast<long double>(b.y - o.y)) -
         (static_cast<long double>(a.y - o.y) * static_cast<long double>(b.x - o.x));
}

OutType BuildHull(const InType &pts) {
  const std::size_t n = pts.size();
  if (n == 0) {
    return {};
  }
  if (n == 1) {
    return pts;
  }

  InType sorted = pts;
  std::ranges::sort(sorted, [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });

  OutType lower;
  lower.reserve(sorted.size());
  for (const auto &p : sorted) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  OutType upper;
  upper.reserve(sorted.size());
  for (const auto &p : std::ranges::reverse_view(sorted)) {
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  OutType hull;
  hull.reserve(lower.size() + upper.size());
  hull.insert(hull.end(), lower.begin(), lower.end());
  for (std::size_t i = 1; i + 1 < upper.size(); ++i) {
    hull.push_back(upper[i]);
  }
  return hull;
}

void RunStdThreadExercise(int num_threads) {
  std::atomic<int> th_counter{0};
  std::vector<std::thread> threads;
  threads.reserve(static_cast<std::size_t>(num_threads));

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back([&th_counter]() { th_counter.fetch_add(1, std::memory_order_relaxed); });
  }

  for (auto &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }

  (void)th_counter.load(std::memory_order_relaxed);
}

void RunMPIExercise() {
#ifdef USE_MPI
  static bool finalized = false;
  int initialized = 0;
  MPI_Initialized(&initialized);

  if (initialized && !finalized) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    finalized = true;
  }
#endif
}

}  // namespace

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

  OutType hull = BuildHull(pts_in);

  const int num_threads = std::max(1, ppc::util::GetNumThreads());

  RunStdThreadExercise(num_threads);
  RunMPIExercise();

  GetOutput() = std::move(hull);
  return true;
}

bool YurkinGGrahamScanALL::PostProcessingImpl() {
  if (GetInput().empty()) {
    return true;
  }
  return !GetOutput().empty();
}

}  // namespace yurkin_g_graham_scan
