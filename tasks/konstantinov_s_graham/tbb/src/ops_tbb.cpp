#include "konstantinov_s_graham/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <utility>
#include <vector>

#include "konstantinov_s_graham/common/include/common.hpp"

namespace konstantinov_s_graham {

bool KonstantinovAGrahamTBB::IsLowerAnchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t lhs,
                                           size_t rhs) {
  if (ys[lhs] < ys[rhs] - kKEps) {
    return true;
  }

  if (std::abs(ys[lhs] - ys[rhs]) < kKEps && xs[lhs] < xs[rhs] - kKEps) {
    return true;
  }

  return false;
}

KonstantinovAGrahamTBB::KonstantinovAGrahamTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KonstantinovAGrahamTBB::ValidationImpl() {
  return GetInput().first.size() == GetInput().second.size();
}

bool KonstantinovAGrahamTBB::PreProcessingImpl() {
  return true;
}

void KonstantinovAGrahamTBB::RemoveDuplicates(std::vector<double> &xs, std::vector<double> &ys) {
  std::vector<std::pair<double, double>> pts;
  pts.reserve(xs.size());

  for (size_t i = 0; i < xs.size(); ++i) {
    pts.emplace_back(xs[i], ys[i]);
  }

  std::ranges::sort(pts, [](const auto &a, const auto &b) {
    if (std::abs(a.first - b.first) > kKEps) {
      return a.first < b.first;
    }
    return a.second < b.second;
  });

  const auto new_end = std::ranges::unique(pts, [](const auto &a, const auto &b) {
    return std::abs(a.first - b.first) < kKEps && std::abs(a.second - b.second) < kKEps;
  });

  pts.erase(new_end.begin(), pts.end());

  xs.resize(pts.size());
  ys.resize(pts.size());

  for (size_t i = 0; i < pts.size(); ++i) {
    xs[i] = pts[i].first;
    ys[i] = pts[i].second;
  }
}

size_t KonstantinovAGrahamTBB::FindAnchorIndex(const std::vector<double> &xs, const std::vector<double> &ys) {
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

double KonstantinovAGrahamTBB::Dist2(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j) {
  const double dx = xs[j] - xs[i];
  const double dy = ys[j] - ys[i];
  return (dx * dx) + (dy * dy);
}

double KonstantinovAGrahamTBB::CrossVal(const std::vector<double> &xs, const std::vector<double> &ys, size_t i,
                                        size_t j, size_t k) {
  const double ax = xs[j] - xs[i];
  const double ay = ys[j] - ys[i];
  const double bx = xs[k] - xs[i];
  const double by = ys[k] - ys[i];
  return (ax * by) - (ay * bx);
}

std::vector<size_t> KonstantinovAGrahamTBB::CollectAndSortIndices(const std::vector<double> &xs,
                                                                  const std::vector<double> &ys, size_t anchor_idx) {
  std::vector<size_t> idxs(xs.size() - 1);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, xs.size()), [&idxs, anchor_idx](const tbb::blocked_range<size_t> &r) {
    for (size_t i = r.begin(); i < r.end(); ++i) {
      if (i < anchor_idx) {
        idxs[i] = i;
      } else if (i > anchor_idx) {
        idxs[i - 1] = i;
      }
    }
  });

  std::ranges::sort(idxs, [&xs, &ys, anchor_idx](size_t a, size_t b) {
    const double cr = CrossVal(xs, ys, anchor_idx, a, b);

    if (std::abs(cr) < kKEps) {
      return Dist2(xs, ys, anchor_idx, a) < Dist2(xs, ys, anchor_idx, b);
    }

    return cr > 0;
  });

  return idxs;
}

bool KonstantinovAGrahamTBB::AllCollinearWithAnchor(const std::vector<double> &xs, const std::vector<double> &ys,
                                                    size_t anchor_idx, const std::vector<size_t> &sorted_idxs) {
  if (sorted_idxs.empty()) {
    return true;
  }

  return tbb::parallel_reduce(tbb::blocked_range<size_t>(1, sorted_idxs.size()), true,
                              [&xs, &ys, anchor_idx, &sorted_idxs](const tbb::blocked_range<size_t> &r, bool local_ok) {
    for (size_t i = r.begin(); i < r.end() && local_ok; ++i) {
      if (std::abs(CrossVal(xs, ys, anchor_idx, sorted_idxs[0], sorted_idxs[i])) > kKEps) {
        local_ok = false;
      }
    }
    return local_ok;
  }, [](bool lhs, bool rhs) { return lhs && rhs; });
}

std::vector<std::pair<double, double>> KonstantinovAGrahamTBB::BuildHullFromSorted(
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

bool KonstantinovAGrahamTBB::RunImpl() {
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

bool KonstantinovAGrahamTBB::PostProcessingImpl() {
  return true;
}

}  // namespace konstantinov_s_graham
