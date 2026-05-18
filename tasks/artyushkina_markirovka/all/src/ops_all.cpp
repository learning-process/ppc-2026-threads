#include "artyushkina_markirovka/all/include/ops_all.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"

namespace artyushkina_markirovka {

struct NeighborOffsetAll {
  int di;
  int dj;
  int check_i_min;
  int check_i_max;
  int check_j_min;
  int check_j_max;
};

namespace {

struct NeighborInfo {
  int ni;
  int nj;
  int label;
};

std::vector<NeighborOffsetAll> GetFirstPassNeighbors() {
  std::vector<NeighborOffsetAll> neighbors(4);
  neighbors[0] = {.di = -1, .dj = -1, .check_i_min = 1, .check_i_max = 0, .check_j_min = 1, .check_j_max = 0};
  neighbors[1] = {.di = -1, .dj = 0, .check_i_min = 1, .check_i_max = 0, .check_j_min = 0, .check_j_max = 0};
  neighbors[2] = {.di = -1, .dj = 1, .check_i_min = 1, .check_i_max = 0, .check_j_min = 0, .check_j_max = 1};
  neighbors[3] = {.di = 0, .dj = -1, .check_i_min = 0, .check_i_max = 0, .check_j_min = 1, .check_j_max = 0};
  return neighbors;
}

bool IsValidNeighborOffset(int i, int j, int rows, int cols, const NeighborOffsetAll &offset) {
  if (offset.check_i_min == 1 && i <= 0) {
    return false;
  }
  if (offset.check_i_max == 1 && i >= rows - 1) {
    return false;
  }
  if (offset.check_j_min == 1 && j <= 0) {
    return false;
  }
  if (offset.check_j_max == 1 && j >= cols - 1) {
    return false;
  }
  return true;
}

NeighborInfo GetNeighborInfo(int i, int j, int rows, int cols, const std::vector<std::vector<int>> &labels,
                             const NeighborOffsetAll &offset) {
  NeighborInfo info = {.ni = -1, .nj = -1, .label = 0};
  if (!IsValidNeighborOffset(i, j, rows, cols, offset)) {
    return info;
  }
  info.ni = i + offset.di;
  info.nj = j + offset.dj;
  info.label = labels[static_cast<size_t>(info.ni)][static_cast<size_t>(info.nj)];
  return info;
}

void FindMinLabelFromNeighbors(int i, int j, int rows, int cols, const std::vector<std::vector<int>> &labels,
                               const std::vector<NeighborOffsetAll> &neighbors, int &min_label, bool &has_neighbors) {
  for (const auto &offset : neighbors) {
    NeighborInfo info = GetNeighborInfo(i, j, rows, cols, labels, offset);
    if (info.label != 0) {
      has_neighbors = true;
      min_label = std::min(info.label, min_label);
    }
  }
}

void UnionNeighborLabels(int i, int j, int rows, int cols, std::vector<std::vector<int>> &labels,
                         std::vector<int> &equivalent_labels, const std::vector<NeighborOffsetAll> &neighbors,
                         int min_label) {
  for (const auto &offset : neighbors) {
    NeighborInfo info = GetNeighborInfo(i, j, rows, cols, labels, offset);
    if (info.label != 0 && info.label != min_label) {
      MarkingComponentsALL::UnionLabels(equivalent_labels, info.label, min_label);
    }
  }
}

void AssignNewLabel(int i, int j, int &next_label, std::vector<std::vector<int>> &labels,
                    std::vector<int> &equivalent_labels) {
  labels[static_cast<size_t>(i)][static_cast<size_t>(j)] = next_label;
  equivalent_labels.push_back(next_label);
  ++next_label;
}

void AssignExistingLabel(int i, int j, int min_label, std::vector<std::vector<int>> &labels) {
  labels[static_cast<size_t>(i)][static_cast<size_t>(j)] = min_label;
}

}  // namespace

MarkingComponentsALL::MarkingComponentsALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MarkingComponentsALL::ValidationImpl() {
  return GetInput().size() >= 2;
}

bool MarkingComponentsALL::PreProcessingImpl() {
  const auto &input = GetInput();
  rows_ = static_cast<int>(input[0]);
  cols_ = static_cast<int>(input[1]);

  labels_.clear();
  labels_.resize(static_cast<size_t>(rows_));
  for (int i = 0; i < rows_; ++i) {
    labels_[static_cast<size_t>(i)].resize(static_cast<size_t>(cols_), 0);
  }

  equivalent_labels_.clear();
  equivalent_labels_.push_back(0);

  return true;
}

int MarkingComponentsALL::FindRoot(std::vector<int> &parent, int label) {
  int root = label;
  while (parent[static_cast<size_t>(root)] != root) {
    root = parent[static_cast<size_t>(root)];
  }
  int current = label;
  while (current != root) {
    int next = parent[static_cast<size_t>(current)];
    parent[static_cast<size_t>(current)] = root;
    current = next;
  }
  return root;
}

void MarkingComponentsALL::UnionLabels(std::vector<int> &parent, int label1, int label2) {
  int root1 = FindRoot(parent, label1);
  int root2 = FindRoot(parent, label2);
  if (root1 != root2) {
    if (root1 < root2) {
      parent[static_cast<size_t>(root2)] = root1;
    } else {
      parent[static_cast<size_t>(root1)] = root2;
    }
  }
}

void MarkingComponentsALL::FirstPass() {
  const auto &input = GetInput();
  int next_label = 1;
  auto neighbors = GetFirstPassNeighbors();

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      size_t idx = (static_cast<size_t>(i) * static_cast<size_t>(cols_)) + static_cast<size_t>(j) + 2;

      // ИСПРАВЛЕНИЕ: обрабатываем объекты (0), пропускаем фон (255)
      if (input[idx] == 255) {
        continue;
      }

      int min_label = next_label;
      bool has_neighbors = false;
      FindMinLabelFromNeighbors(i, j, rows_, cols_, labels_, neighbors, min_label, has_neighbors);

      if (!has_neighbors) {
        AssignNewLabel(i, j, next_label, labels_, equivalent_labels_);
      } else {
        AssignExistingLabel(i, j, min_label, labels_);
        UnionNeighborLabels(i, j, rows_, cols_, labels_, equivalent_labels_, neighbors, min_label);
      }
    }
  }
}

void MarkingComponentsALL::SecondPass() {
  int label_count = static_cast<int>(equivalent_labels_.size());

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (labels_[static_cast<size_t>(i)][static_cast<size_t>(j)] != 0) {
        labels_[static_cast<size_t>(i)][static_cast<size_t>(j)] =
            FindRoot(equivalent_labels_, labels_[static_cast<size_t>(i)][static_cast<size_t>(j)]);
      }
    }
  }

  std::vector<int> remap(static_cast<size_t>(label_count), 0);
  int current_label = 1;
  for (int i = 1; i < label_count; ++i) {
    if (equivalent_labels_[static_cast<size_t>(i)] == i) {
      remap[static_cast<size_t>(i)] = current_label;
      ++current_label;
    }
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (labels_[static_cast<size_t>(i)][static_cast<size_t>(j)] != 0) {
        labels_[static_cast<size_t>(i)][static_cast<size_t>(j)] =
            remap[static_cast<size_t>(labels_[static_cast<size_t>(i)][static_cast<size_t>(j)])];
      }
    }
  }
}

bool MarkingComponentsALL::RunImpl() {
  const auto &input = GetInput();
  if (input.size() < 2 || rows_ == 0 || cols_ == 0) {
    return false;
  }
  FirstPass();
  SecondPass();
  return true;
}

bool MarkingComponentsALL::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();
  output.reserve((static_cast<size_t>(rows_) * static_cast<size_t>(cols_)) + 2);
  output.push_back(rows_);
  output.push_back(cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      output.push_back(labels_[static_cast<size_t>(i)][static_cast<size_t>(j)]);
    }
  }
  return true;
}

}  // namespace artyushkina_markirovka
