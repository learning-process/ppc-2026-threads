#include "konstantinov_s_graham/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <thread>
#include <utility>
#include <vector>

#include "konstantinov_s_graham/common/include/common.hpp"

namespace konstantinov_s_graham {

bool KonstantinovAGrahamSTL::IsLowerAnchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t lhs,
                                           size_t rhs) {
  if (ys[lhs] < ys[rhs] - kKEps) {
    return true;
  }

  if (std::abs(ys[lhs] - ys[rhs]) < kKEps && xs[lhs] < xs[rhs] - kKEps) {
    return true;
  }

  return false;
}

size_t KonstantinovAGrahamSTL::GetThreadCount(size_t n) {
  const size_t hw = std::max<size_t>(1, static_cast<size_t>(std::thread::hardware_concurrency()));
  return std::min(hw, n);
}

size_t KonstantinovAGrahamSTL::FindLocalAnchor(const std::vector<double> &xs, const std::vector<double> &ys,
                                               size_t begin, size_t end) {
  size_t best = begin;

  for (size_t i = begin + 1; i < end; ++i) {
    if (IsLowerAnchor(xs, ys, i, best)) {
      best = i;
    }
  }

  return best;
}

KonstantinovAGrahamSTL::KonstantinovAGrahamSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KonstantinovAGrahamSTL::ValidationImpl() {
  return GetInput().first.size() == GetInput().second.size();
}

bool KonstantinovAGrahamSTL::PreProcessingImpl() {
  return true;
}

void KonstantinovAGrahamSTL::RemoveDuplicates(std::vector<double> &xs, std::vector<double> &ys) {
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

size_t KonstantinovAGrahamSTL::FindAnchorIndex(const std::vector<double> &xs, const std::vector<double> &ys) {
  if (xs.empty()) {
    return 0;
  }

  const size_t thread_count = GetThreadCount(xs.size());

  if (thread_count <= 1 || xs.size() < (thread_count * 2)) {
    return FindLocalAnchor(xs, ys, 0, xs.size());
  }

  std::vector<size_t> local_results(thread_count);
  std::vector<std::thread> workers;
  workers.reserve(thread_count);

  const size_t block_size = (xs.size() + thread_count - 1) / thread_count;

  for (size_t thread_idx = 0; thread_idx < thread_count; ++thread_idx) {
    const size_t begin = thread_idx * block_size;
    const size_t end = std::min(xs.size(), begin + block_size);

    if (begin >= end) {
      local_results[thread_idx] = 0;
      continue;
    }

    workers.emplace_back(
        [&, thread_idx, begin, end]() { local_results[thread_idx] = FindLocalAnchor(xs, ys, begin, end); });
  }

  for (auto &worker : workers) {
    worker.join();
  }

  size_t best = local_results[0];

  for (size_t thread_idx = 1; thread_idx < thread_count; ++thread_idx) {
    if (IsLowerAnchor(xs, ys, local_results[thread_idx], best)) {
      best = local_results[thread_idx];
    }
  }

  return best;
}

double KonstantinovAGrahamSTL::Dist2(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j) {
  const double dx = xs[j] - xs[i];
  const double dy = ys[j] - ys[i];
  return (dx * dx) + (dy * dy);
}

double KonstantinovAGrahamSTL::CrossVal(const std::vector<double> &xs, const std::vector<double> &ys, size_t i,
                                        size_t j, size_t k) {
  const double ax = xs[j] - xs[i];
  const double ay = ys[j] - ys[i];
  const double bx = xs[k] - xs[i];
  const double by = ys[k] - ys[i];
  return (ax * by) - (ay * bx);
}

std::vector<size_t> KonstantinovAGrahamSTL::CollectAndSortIndices(const std::vector<double> &xs,
                                                                  const std::vector<double> &ys, size_t anchor_idx) {
  std::vector<size_t> idxs(xs.size() - 1);

  const size_t thread_count = GetThreadCount(xs.size());

  if (thread_count <= 1 || xs.size() < (thread_count * 2)) {
    FillIndexRange(idxs, 0, xs.size(), anchor_idx);
  } else {
    FillIndicesParallel(idxs, xs.size(), anchor_idx, thread_count);
  }

  std::ranges::sort(idxs, [&xs, &ys, anchor_idx](size_t lhs, size_t rhs) {
    const double cross = CrossVal(xs, ys, anchor_idx, lhs, rhs);

    if (std::abs(cross) < kKEps) {
      return Dist2(xs, ys, anchor_idx, lhs) < Dist2(xs, ys, anchor_idx, rhs);
    }

    return cross > 0;
  });

  return idxs;
}

bool KonstantinovAGrahamSTL::CheckCollinearRange(const std::vector<double> &xs, const std::vector<double> &ys,
                                                 size_t anchor_idx, const std::vector<size_t> &sorted_idxs,
                                                 size_t begin, size_t end) {
  for (size_t i = begin; i < end; ++i) {
    if (std::abs(CrossVal(xs, ys, anchor_idx, sorted_idxs[0], sorted_idxs[i])) > kKEps) {
      return false;
    }
  }

  return true;
}

void KonstantinovAGrahamSTL::FillIndexRange(std::vector<size_t> &idxs, size_t begin, size_t end, size_t anchor_idx) {
  for (size_t i = begin; i < end; ++i) {
    if (i < anchor_idx) {
      idxs[i] = i;
    } else if (i > anchor_idx) {
      idxs[i - 1] = i;
    }
  }
}

void KonstantinovAGrahamSTL::FillIndicesParallel(std::vector<size_t> &idxs, size_t point_count, size_t anchor_idx,
                                                 size_t thread_count) {
  std::vector<std::thread> workers;
  workers.reserve(thread_count);

  const size_t block_size = (point_count + thread_count - 1) / thread_count;

  for (size_t thread_idx = 0; thread_idx < thread_count; ++thread_idx) {
    const size_t begin = thread_idx * block_size;
    const size_t end = std::min(point_count, begin + block_size);

    if (begin >= end) {
      continue;
    }

    workers.emplace_back([&, begin, end]() { FillIndexRange(idxs, begin, end, anchor_idx); });
  }

  for (auto &worker : workers) {
    worker.join();
  }
}

bool KonstantinovAGrahamSTL::AllCollinearWithAnchor(const std::vector<double> &xs, const std::vector<double> &ys,
                                                    size_t anchor_idx, const std::vector<size_t> &sorted_idxs) {
  if (sorted_idxs.empty()) {
    return true;
  }

  const size_t thread_count = GetThreadCount(sorted_idxs.size());

  if (thread_count <= 1 || sorted_idxs.size() < (thread_count * 2)) {
    return CheckCollinearRange(xs, ys, anchor_idx, sorted_idxs, 1, sorted_idxs.size());
  }

  std::atomic<bool> result{true};

  std::vector<std::thread> workers;
  workers.reserve(thread_count);

  const size_t block_size = (sorted_idxs.size() + thread_count - 1) / thread_count;

  for (size_t thread_idx = 0; thread_idx < thread_count; ++thread_idx) {
    const size_t begin = std::max<size_t>(1, thread_idx * block_size);
    const size_t end = std::min(sorted_idxs.size(), begin + block_size);

    if (begin >= end) {
      continue;
    }

    workers.emplace_back([&, begin, end]() {
      if (!CheckCollinearRange(xs, ys, anchor_idx, sorted_idxs, begin, end)) {
        result.store(false);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }

  return result.load();
}

std::vector<std::pair<double, double>> KonstantinovAGrahamSTL::BuildHullFromSorted(
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
      const double cr = CrossVal(xs, ys, p, q, cur);

      if (cr <= kKEps) {
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

bool KonstantinovAGrahamSTL::RunImpl() {
  const InType &inp = GetInput();
  auto xs = inp.first;
  auto ys = inp.second;

  RemoveDuplicates(xs, ys);

  if (xs.size() != ys.size() || xs.empty()) {
    GetOutput() = {};
    return true;
  }

  if (xs.size() < 3) {
    std::vector<std::pair<double, double>> out;
    out.reserve(xs.size());

    for (size_t i = 0; i < xs.size(); ++i) {
      out.emplace_back(xs[i], ys[i]);
    }

    GetOutput() = out;
    return true;
  }

  const size_t anchor = FindAnchorIndex(xs, ys);
  std::vector<size_t> sorted_idxs = CollectAndSortIndices(xs, ys, anchor);

  if (sorted_idxs.empty()) {
    GetOutput() = {{xs[anchor], ys[anchor]}};
    return true;
  }

  if (AllCollinearWithAnchor(xs, ys, anchor, sorted_idxs)) {
    const size_t far_idx = sorted_idxs.back();
    GetOutput() = {{xs[anchor], ys[anchor]}, {xs[far_idx], ys[far_idx]}};
    return true;
  }

  GetOutput() = BuildHullFromSorted(xs, ys, anchor, sorted_idxs);
  return true;
}

bool KonstantinovAGrahamSTL::PostProcessingImpl() {
  return true;
}

}  // namespace konstantinov_s_graham
