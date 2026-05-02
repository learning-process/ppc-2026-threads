#include "dergachev_a_graham_scan/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <thread>
#include <utility>
#include <vector>

#include "dergachev_a_graham_scan/common/include/common.hpp"
#include "util/include/util.hpp"

namespace dergachev_a_graham_scan {

namespace {

using Pt = std::pair<double, double>;

double Cross(const Pt &o, const Pt &a, const Pt &b) {
  return ((a.first - o.first) * (b.second - o.second)) - ((a.second - o.second) * (b.first - o.first));
}

double Dist2(const Pt &a, const Pt &b) {
  double dx = a.first - b.first;
  double dy = a.second - b.second;
  return (dx * dx) + (dy * dy);
}

const double kPi = std::acos(-1.0);

int FindPivot(const std::vector<Pt> &pts) {
  int best = 0;
  for (int i = 1; std::cmp_less(i, pts.size()); i++) {
    if (pts[i].second < pts[best].second || (pts[i].second == pts[best].second && pts[i].first < pts[best].first)) {
      best = i;
    }
  }
  return best;
}

bool AngleCompare(const Pt &pivot, const Pt &a, const Pt &b) {
  double c = Cross(pivot, a, b);
  if (c != 0.0) {
    return c > 0.0;
  }
  return Dist2(pivot, a) < Dist2(pivot, b);
}

void AngleSort(std::vector<Pt> &pts) {
  Pt pivot = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [&pivot](const Pt &a, const Pt &b) { return AngleCompare(pivot, a, b); });
}

std::vector<Pt> BuildHull(std::vector<Pt> pts) {
  if (pts.size() < 2) {
    return pts;
  }
  int pivot = FindPivot(pts);
  std::swap(pts[0], pts[pivot]);
  AngleSort(pts);
  std::vector<Pt> hull;
  for (const auto &p : pts) {
    while (hull.size() > 1 && Cross(hull[hull.size() - 2], hull.back(), p) <= 0.0) {
      hull.pop_back();
    }
    hull.push_back(p);
  }
  return hull;
}

int ChunkLen(int idx, int total, int parts) {
  return (total / parts) + ((idx < (total % parts)) ? 1 : 0);
}

std::vector<Pt> ThreadedHull(const std::vector<Pt> &pts) {
  int n = static_cast<int>(pts.size());
  int num_threads = ppc::util::GetNumThreads();
  if (n < num_threads * 4) {
    return BuildHull({pts.begin(), pts.end()});
  }
  std::vector<std::vector<Pt>> partial(num_threads);
  std::vector<std::thread> workers;
  int off = 0;
  for (int ti = 0; ti < num_threads; ti++) {
    int len = ChunkLen(ti, n, num_threads);
    workers.emplace_back(
        [&partial, &pts, off, len, ti]() { partial[ti] = BuildHull({pts.begin() + off, pts.begin() + off + len}); });
    off += len;
  }
  for (auto &w : workers) {
    w.join();
  }
  std::vector<Pt> merged;
  for (const auto &h : partial) {
    merged.insert(merged.end(), h.begin(), h.end());
  }
  return BuildHull(std::move(merged));
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
  int n = GetInput();
  if (n <= 0) {
    points_.clear();
    return true;
  }
  points_.resize(n);
  double step = (2.0 * kPi) / n;
  for (int i = 0; i < n; i++) {
    points_[i] = {std::cos(step * i), std::sin(step * i)};
  }
  if (n > 3) {
    points_.emplace_back(0.0, 0.0);
  }
  return true;
}

bool DergachevAGrahamScanALL::RunImpl() {
  hull_.clear();
  std::vector<Pt> pts(points_.begin(), points_.end());

  if (pts.size() <= 1) {
    hull_ = pts;
    return true;
  }

  hull_ = ThreadedHull(pts);
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool DergachevAGrahamScanALL::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan
