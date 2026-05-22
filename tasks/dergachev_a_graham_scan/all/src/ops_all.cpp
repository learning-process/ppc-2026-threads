#include "dergachev_a_graham_scan/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
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

void AngleSort(std::vector<Pt> &pts) {
  Pt pivot = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [&pivot](const Pt &a, const Pt &b) {
    double c = Cross(pivot, a, b);
    if (c > 0.0) {
      return true;
    }
    if (c < 0.0) {
      return false;
    }
    return Dist2(pivot, a) < Dist2(pivot, b);
  });
}

std::vector<Pt> BuildHull(std::vector<Pt> pts) {
  int n = static_cast<int>(pts.size());
  if (n <= 1) {
    return pts;
  }
  if (std::all_of(pts.begin() + 1, pts.end(),
                  [&](const Pt &p) { return p.first == pts[0].first && p.second == pts[0].second; })) {
    return {pts[0]};
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

std::vector<Pt> ThreadedHull(const std::vector<Pt> &pts, int num_threads) {
  int n = static_cast<int>(pts.size());
  if (num_threads <= 1 || n < num_threads * 4) {
    return BuildHull({pts.begin(), pts.end()});
  }
  std::vector<std::vector<Pt>> partial(num_threads);
  std::vector<std::thread> workers;
  int off = 0;
  for (int ti = 0; ti < num_threads; ti++) {
    int len = (n / num_threads) + ((ti < (n % num_threads)) ? 1 : 0);
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

std::vector<double> Flatten(const std::vector<Pt> &pts) {
  std::vector<double> flat(pts.size() * 2);
  for (size_t i = 0; i < pts.size(); i++) {
    flat[i * 2] = pts[i].first;
    flat[(i * 2) + 1] = pts[i].second;
  }
  return flat;
}

std::vector<Pt> Unflatten(const std::vector<double> &flat) {
  std::vector<Pt> pts(flat.size() / 2);
  for (size_t i = 0; i < pts.size(); i++) {
    pts[i] = {flat[i * 2], flat[(i * 2) + 1]};
  }
  return pts;
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
  int n = static_cast<int>(points_.size());

  if (n <= 1) {
    if (!points_.empty()) {
      hull_.push_back(points_[0]);
    }
    return true;
  }

  if (std::all_of(points_.begin() + 1, points_.end(),
                  [&](const Pt &p) { return p.first == points_[0].first && p.second == points_[0].second; })) {
    hull_.push_back(points_[0]);
    return true;
  }

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> send_counts(world_size);
  std::vector<int> send_displs(world_size);
  int disp = 0;
  for (int i = 0; i < world_size; i++) {
    int chunk = (n / world_size) + ((i < (n % world_size)) ? 1 : 0);
    send_counts[i] = chunk * 2;
    send_displs[i] = disp;
    disp += send_counts[i];
  }

  std::vector<double> flat_all;
  if (rank == 0) {
    flat_all = Flatten(points_);
  }

  int local_size = send_counts[rank];
  std::vector<double> local_flat(local_size);
  MPI_Scatterv(flat_all.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE, local_flat.data(), local_size,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::vector<Pt> local_pts = Unflatten(local_flat);
  int num_threads = ppc::util::GetNumThreads();
  std::vector<Pt> local_hull = ThreadedHull(local_pts, num_threads);

  int local_hull_flat_size = static_cast<int>(local_hull.size()) * 2;
  std::vector<int> recv_counts(world_size);
  MPI_Gather(&local_hull_flat_size, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> recv_displs(world_size);
  int total_recv = 0;
  if (rank == 0) {
    for (int i = 0; i < world_size; i++) {
      recv_displs[i] = total_recv;
      total_recv += recv_counts[i];
    }
  }

  std::vector<double> hull_flat = Flatten(local_hull);
  std::vector<double> gathered_flat(total_recv);
  MPI_Gatherv(hull_flat.data(), local_hull_flat_size, MPI_DOUBLE, gathered_flat.data(), recv_counts.data(),
              recv_displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::vector<Pt> all_hull_pts = Unflatten(gathered_flat);
    hull_ = BuildHull(std::move(all_hull_pts));
  }

  int hull_size = static_cast<int>(hull_.size());
  MPI_Bcast(&hull_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  hull_.resize(hull_size);
  std::vector<double> result_flat(static_cast<size_t>(hull_size) * 2);
  if (rank == 0) {
    result_flat = Flatten(hull_);
  }
  MPI_Bcast(result_flat.data(), hull_size * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    hull_ = Unflatten(result_flat);
  }

  return true;
}

bool DergachevAGrahamScanALL::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan
