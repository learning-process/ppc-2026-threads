#include "makovskiy_i_graham_hull/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <future>
#include <utility>
#include <vector>

#include "makovskiy_i_graham_hull/common/include/common.hpp"
#include "util/include/util.hpp"

namespace makovskiy_i_graham_hull {

namespace {

double CrossProduct(const Point &o, const Point &a, const Point &b) {
  return ((a.x - o.x) * (b.y - o.y)) - ((a.y - o.y) * (b.x - o.x));
}

double DistSq(const Point &a, const Point &b) {
  return ((a.x - b.x) * (a.x - b.x)) + ((a.y - b.y) * (a.y - b.y));
}

bool IsBetterMin(const Point &candidate, const Point &current_min) {
  if (candidate.y < current_min.y - 1e-9) {
    return true;
  }
  return (std::abs(candidate.y - current_min.y) <= 1e-9) && (candidate.x < current_min.x);
}

size_t FindMinPointIndexSTL(const std::vector<Point> &points) {
  const size_t n = points.size();
  const int num_threads = std::max(1, ppc::util::GetNumThreads());

  std::vector<std::future<size_t>> futures;
  const size_t chunk = (n + static_cast<size_t>(num_threads) - 1) / static_cast<size_t>(num_threads);

  auto worker = [&points](size_t start, size_t end) {
    size_t local_min = start;
    for (size_t j = start + 1; j < end; ++j) {
      if (IsBetterMin(points[j], points[local_min])) {
        local_min = j;
      }
    }
    return local_min;
  };

  for (int i = 0; i < num_threads; ++i) {
    const size_t start = static_cast<size_t>(i) * chunk;
    const size_t end = std::min(start + chunk, n);
    if (start >= n) {
      break;
    }
    futures.push_back(std::async(std::launch::async, worker, start, end));
  }

  size_t min_idx = futures[0].get();
  for (size_t i = 1; i < futures.size(); ++i) {
    const size_t local_min = futures[i].get();
    if (IsBetterMin(points[local_min], points[min_idx])) {
      min_idx = local_min;
    }
  }

  return min_idx;
}

template <typename RandomIt, typename Compare>
void StlParallelSortSub(RandomIt first, RandomIt last, Compare comp) {
  if (last - first < 32) {
    std::sort(first, last, comp);
    return;
  }
  const auto pivot = *(first + ((last - first) / 2));
  RandomIt middle1 = std::partition(first, last, [pivot, comp](const auto &a) { return comp(a, pivot); });
  RandomIt middle2 = std::partition(middle1, last, [pivot, comp](const auto &a) { return !comp(pivot, a); });

  auto future1 = std::async(std::launch::async, [first, middle1, comp]() { std::sort(first, middle1, comp); });
  std::sort(middle2, last, comp);
  future1.wait();
}

template <typename RandomIt, typename Compare>
void StlParallelSort(RandomIt first, RandomIt last, Compare comp) {
  if (last - first < 32) {
    std::sort(first, last, comp);
    return;
  }
  const auto pivot = *(first + ((last - first) / 2));
  RandomIt middle1 = std::partition(first, last, [pivot, comp](const auto &a) { return comp(a, pivot); });
  RandomIt middle2 = std::partition(middle1, last, [pivot, comp](const auto &a) { return !comp(pivot, a); });

  auto future1 = std::async(std::launch::async, [first, middle1, comp]() { StlParallelSortSub(first, middle1, comp); });
  StlParallelSortSub(middle2, last, comp);
  future1.wait();
}

std::vector<Point> FilterPointsSTL(const std::vector<Point> &points, const Point &p0) {
  const size_t n = points.size();
  if (n <= 2) {
    return points;
  }

  std::vector<uint8_t> keep(n, 1);
  const int num_threads = std::max(1, ppc::util::GetNumThreads());

  const size_t num_elements = n - 2;
  const size_t chunk = (num_elements + static_cast<size_t>(num_threads) - 1) / static_cast<size_t>(num_threads);
  std::vector<std::future<void>> futures;

  auto worker = [&points, &keep, p0](size_t start, size_t end) {
    for (size_t j = start; j < end; ++j) {
      if (std::abs(CrossProduct(p0, points[j], points[j + 1])) < 1e-9) {
        keep[j] = 0;
      }
    }
  };

  for (int i = 0; i < num_threads; ++i) {
    const size_t start = 1 + (static_cast<size_t>(i) * chunk);
    const size_t end = std::min(start + chunk, n - 1);

    if (start >= n - 1) {
      break;
    }

    futures.push_back(std::async(std::launch::async, worker, start, end));
  }

  for (auto &f : futures) {
    f.wait();
  }

  std::vector<Point> filtered;
  filtered.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (keep[i] != 0) {
      filtered.push_back(points[i]);
    }
  }
  return filtered;
}

std::vector<Point> BuildHull(const std::vector<Point> &filtered) {
  std::vector<Point> hull;
  hull.push_back(filtered[0]);
  hull.push_back(filtered[1]);
  hull.push_back(filtered[2]);

  for (size_t i = 3; i < filtered.size(); ++i) {
    while (hull.size() > 1 && CrossProduct(hull[hull.size() - 2], hull.back(), filtered[i]) <= 1e-9) {
      hull.pop_back();
    }
    hull.push_back(filtered[i]);
  }
  return hull;
}

std::vector<Point> ComputeHullSTL(const std::vector<Point> &input_points) {
  if (input_points.size() < 3) {
    return input_points;
  }

  std::vector<Point> points = input_points;
  const size_t min_idx = FindMinPointIndexSTL(points);

  std::swap(points[0], points[min_idx]);
  const Point p0 = points[0];

  auto comp = [p0](const Point &a, const Point &b) {
    const double cp = CrossProduct(p0, a, b);
    if (std::abs(cp) < 1e-9) {
      return DistSq(p0, a) < DistSq(p0, b);
    }
    return cp > 0;
  };

  StlParallelSort(points.begin() + 1, points.end(), comp);

  std::vector<Point> filtered = FilterPointsSTL(points, p0);

  if (filtered.size() < 3) {
    return filtered;
  }

  return BuildHull(filtered);
}

std::vector<Point> ReceivePointsWorker() {
  std::vector<Point> local_points;
  int count = 0;
  MPI_Recv(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  if (count > 0) {
    local_points.resize(static_cast<size_t>(count));
    MPI_Recv(local_points.data(), count * 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  return local_points;
}

std::vector<Point> SendPointsRoot(int size, const std::vector<Point> &input_points) {
  const int num_points = static_cast<int>(input_points.size());

  const int chunk = num_points / size;
  const int rem = num_points % size;

  int offset = chunk + (rem > 0 ? 1 : 0);
  for (int i = 1; i < size; ++i) {
    const int count = chunk + (i < rem ? 1 : 0);
    MPI_Send(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    if (count > 0) {
      MPI_Send(input_points.data() + offset, count * 2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
    }
    offset += count;
  }

  const int local_count = chunk + (rem > 0 ? 1 : 0);
  std::vector<Point> local_points;
  local_points.assign(input_points.begin(), input_points.begin() + local_count);
  return local_points;
}

std::vector<Point> ScatterPoints(int rank, int size, const std::vector<Point> &input_points) {
  if (rank == 0) {
    return SendPointsRoot(size, input_points);
  }
  return ReceivePointsWorker();
}

std::vector<Point> GatherHulls(int rank, int size, const std::vector<Point> &local_hull) {
  int local_hull_size = static_cast<int>(local_hull.size());
  std::vector<int> hull_sizes(static_cast<size_t>(size), 0);

  MPI_Gather(&local_hull_size, 1, MPI_INT, rank == 0 ? hull_sizes.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displs;
  std::vector<int> recvcounts;
  int total_hull_points = 0;

  if (rank == 0) {
    displs.resize(static_cast<size_t>(size), 0);
    recvcounts.resize(static_cast<size_t>(size), 0);
    for (int i = 0; i < size; ++i) {
      const auto idx = static_cast<size_t>(i);
      recvcounts[idx] = hull_sizes[idx] * 2;
      displs[idx] = total_hull_points * 2;
      total_hull_points += hull_sizes[idx];
    }
  }

  std::vector<Point> all_hulls;
  if (rank == 0) {
    all_hulls.resize(static_cast<size_t>(total_hull_points));
  }

  MPI_Gatherv(local_hull.data(), local_hull_size * 2, MPI_DOUBLE, rank == 0 ? all_hulls.data() : nullptr,
              recvcounts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return all_hulls;
}

}  // namespace

ConvexHullGrahamALL::ConvexHullGrahamALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ConvexHullGrahamALL::ValidationImpl() {
  return true;
}

bool ConvexHullGrahamALL::PreProcessingImpl() {
  return true;
}

bool ConvexHullGrahamALL::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<Point> input_points;
  if (rank == 0) {
    input_points = GetInput();
  }

  std::vector<Point> local_points = ScatterPoints(rank, size, input_points);
  std::vector<Point> local_hull = ComputeHullSTL(local_points);
  std::vector<Point> all_hulls = GatherHulls(rank, size, local_hull);

  if (rank == 0) {
    GetOutput() = ComputeHullSTL(all_hulls);
  }

  return true;
}

bool ConvexHullGrahamALL::PostProcessingImpl() {
  return true;
}

}  // namespace makovskiy_i_graham_hull
