#include "konstantinov_s_graham/seq/include/ops_seq.hpp"

// #include <numeric>
#include <algorithm>
#include <cstddef>
#include <ranges>
#include <utility>
#include <vector>

#include "konstantinov_s_graham/common/include/common.hpp"
// #include "util/include/util.hpp"

namespace konstantinov_s_graham {

KonstantinovAGrahamSEQ::KonstantinovAGrahamSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  // GetOutput() = 0;
}

bool KonstantinovAGrahamSEQ::ValidationImpl() {
  return GetInput().first.size() == GetInput().second.size();
}

bool KonstantinovAGrahamSEQ::PreProcessingImpl() {
  return true;
}

void KonstantinovAGrahamSEQ::RemoveDuplicates(std::vector<double> &xs, std::vector<double> &ys) {
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

  auto new_end = std::ranges::unique(pts, [](const auto &a, const auto &b) {
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

// Возвращает индекс опорной (самой нижней, при равенстве — самой левой) точки.
size_t KonstantinovAGrahamSEQ::FindAnchorIndex(const std::vector<double> &xs, const std::vector<double> &ys) {
  size_t idx = 0;
  for (size_t i = 1; i < xs.size(); ++i) {
    if (ys[i] < ys[idx] - kKEps || (std::abs(ys[i] - ys[idx]) < kKEps && xs[i] < xs[idx] - kKEps)) {
      idx = i;
    }
  }
  return idx;
}

// Квадрат евклидова расстояния между точками по индексам.
double KonstantinovAGrahamSEQ::Dist2(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j) {
  const double dx = xs[j] - xs[i];
  const double dy = ys[j] - ys[i];
  return (dx * dx) + (dy * dy);
}

// Ориентация троицы точек: >0 левый поворот, <0 правый, 0 коллинеарность.
double KonstantinovAGrahamSEQ::CrossVal(const std::vector<double> &xs, const std::vector<double> &ys, size_t i,
                                        size_t j, size_t k) {
  const double ax = xs[j] - xs[i];
  const double ay = ys[j] - ys[i];
  const double bx = xs[k] - xs[i];
  const double by = ys[k] - ys[i];
  return ax * by - ay * bx;
}

// Подготовить массив индексов всех точек, отсортированных по полярному углу вокруг опорной точки.
// Если углы совпадают — ближняя точка первая.
std::vector<size_t> KonstantinovAGrahamSEQ::CollectAndSortIndices(const std::vector<double> &xs,
                                                                  const std::vector<double> &ys, size_t anchor_idx) {
  std::vector<size_t> idxs;
  idxs.reserve(xs.size() - 1);
  for (size_t i = 0; i < xs.size(); ++i) {
    if (i != anchor_idx) {
      idxs.push_back(i);
    }
  }

  std::ranges::sort(idxs, [&](size_t a, size_t b) {
    double cr = CrossVal(xs, ys, anchor_idx, a, b);

    if (std::abs(cr) < kKEps) {
      return Dist2(xs, ys, anchor_idx, a) < Dist2(xs, ys, anchor_idx, b);
    }

    return cr > 0;
  });

  return idxs;
}

// Проверка, все ли точки (в наборе индексов) коллинеарны с опорной.
bool KonstantinovAGrahamSEQ::AllCollinearWithAnchor(const std::vector<double> &xs, const std::vector<double> &ys,
                                                    size_t anchor_idx, const std::vector<size_t> &sorted_idxs) {
  if (sorted_idxs.empty()) {
    return true;
  }
  for (size_t i = 1; i < sorted_idxs.size(); ++i) {
    if (std::abs(CrossVal(xs, ys, anchor_idx, sorted_idxs[0], sorted_idxs[i])) > kKEps) {
      return false;
    }
  }
  return true;
}

// Построить выпуклую оболочку (возвращает вершины в парах (x,y) по порядку).
std::vector<std::pair<double, double>> KonstantinovAGrahamSEQ::BuildHullFromSorted(
    const std::vector<double> &xs, const std::vector<double> &ys, size_t anchor_idx,
    const std::vector<size_t> &sorted_idxs) {
  std::vector<size_t> stack;
  stack.reserve(sorted_idxs.size() + 1);
  stack.push_back(anchor_idx);
  if (!sorted_idxs.empty()) {
    stack.push_back(sorted_idxs[0]);
  }

  for (size_t i = 1; i < sorted_idxs.size(); ++i) {
    size_t cur = sorted_idxs[i];
    while (stack.size() >= 2) {
      size_t q = stack.back();
      size_t p = stack[stack.size() - 2];
      double cr = CrossVal(xs, ys, p, q, cur);
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

// Вход: два вектора координат одинаковой длины (xs и ys).
// Выход: вектор пар (x,y) — вершины выпуклой оболочки по часовой (или против — зависит от сортировки).
// std::vector<std::pair<double,double>> sequential_graham(const std::vector<double> &xs,
//                                                         const std::vector<double> &ys) {

// }

bool KonstantinovAGrahamSEQ::RunImpl() {
  // std::cout<<"START\n";

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

  size_t anchor = FindAnchorIndex(xs, ys);
  std::vector<size_t> sorted_idxs = CollectAndSortIndices(xs, ys, anchor);
  if (sorted_idxs.empty()) {
    GetOutput() = {{xs[anchor], ys[anchor]}};
    return true;
  }

  if (AllCollinearWithAnchor(xs, ys, anchor, sorted_idxs)) {
    // Возвращаем опорную и самую удалённую вдоль прямой
    size_t far_idx = sorted_idxs.back();
    GetOutput() = {{xs[anchor], ys[anchor]}, {xs[far_idx], ys[far_idx]}};
    return true;
  }

  GetOutput() = BuildHullFromSorted(xs, ys, anchor, sorted_idxs);

  return true;
}

bool KonstantinovAGrahamSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace konstantinov_s_graham
