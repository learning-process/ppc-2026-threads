#include "artyushkina_markirovka/all/include/ops_all.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace artyushkina_markirovka {

// Упрощенная структура без designated initializers
struct NeighborOffsetAll {
  int di;
  int dj;
  int check_i_min;
  int check_i_max;
  int check_j_min;
  int check_j_max;
};

// Инициализация через конструктор вместо designated initializers
static std::vector<NeighborOffsetAll> GetFirstPassNeighbors() {
  std::vector<NeighborOffsetAll> neighbors(4);
  neighbors[0] = {-1, -1, 1, 0, 1, 0};
  neighbors[1] = {-1, 0, 1, 0, 0, 0};
  neighbors[2] = {-1, 1, 1, 0, 0, 1};
  neighbors[3] = {0, -1, 0, 0, 1, 0};
  return neighbors;
}

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
  labels_.resize(rows_);
  for (int i = 0; i < rows_; ++i) {
    labels_[i].resize(cols_, 0);
  }

  equivalent_labels_.clear();
  equivalent_labels_.push_back(0);

  return true;
}

int MarkingComponentsALL::FindRoot(std::vector<int> &parent, int label) {
  while (parent[label] != label) {
    parent[label] = parent[parent[label]];
    label = parent[label];
  }
  return label;
}

void MarkingComponentsALL::UnionLabels(std::vector<int> &parent, int label1, int label2) {
  int root1 = FindRoot(parent, label1);
  int root2 = FindRoot(parent, label2);

  if (root1 != root2) {
    if (root1 < root2) {
      parent[root2] = root1;
    } else {
      parent[root1] = root2;
    }
  }
}

bool MarkingComponentsALL::IsValidNeighbor(int i, int j, const NeighborOffsetAll &offset) const {
  if (offset.check_i_min == 1 && i <= 0) {
    return false;
  }
  if (offset.check_i_max == 1 && i >= rows_ - 1) {
    return false;
  }
  if (offset.check_j_min == 1 && j <= 0) {
    return false;
  }
  if (offset.check_j_max == 1 && j >= cols_ - 1) {
    return false;
  }
  return true;
}

void MarkingComponentsALL::ProcessNeighborFirstPass(int i, int j, const NeighborOffsetAll &offset,
                                                    std::vector<int> &neighbor_labels, int &min_label) {
  if (!IsValidNeighbor(i, j, offset)) {
    return;
  }

  int ni = i + offset.di;
  int nj = j + offset.dj;

  if (labels_[ni][nj] != 0) {
    int label = labels_[ni][nj];
    neighbor_labels.push_back(label);
    if (label < min_label) {
      min_label = label;
    }
  }
}

void MarkingComponentsALL::FirstPass() {
  const auto &input = GetInput();
  int next_label = 1;
  std::vector<NeighborOffsetAll> first_pass_neighbors = GetFirstPassNeighbors();

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      size_t idx = (static_cast<size_t>(i) * static_cast<size_t>(cols_)) + static_cast<size_t>(j) + 2;

      if (input[idx] == 0) {
        std::vector<int> neighbor_labels;
        int min_label = next_label;

        for (size_t k = 0; k < first_pass_neighbors.size(); ++k) {
          ProcessNeighborFirstPass(i, j, first_pass_neighbors[k], neighbor_labels, min_label);
        }

        if (neighbor_labels.empty()) {
          labels_[i][j] = next_label;
          equivalent_labels_.push_back(next_label);
          next_label++;
        } else {
          labels_[i][j] = min_label;
          for (size_t k = 0; k < neighbor_labels.size(); ++k) {
            int label = neighbor_labels[k];
            if (label != min_label) {
              UnionLabels(equivalent_labels_, label, min_label);
            }
          }
        }
      }
    }
  }
}

void MarkingComponentsALL::SecondPass() {
  int label_count = static_cast<int>(equivalent_labels_.size());

  // Первый проход - поиск корней
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (labels_[i][j] != 0) {
        labels_[i][j] = FindRoot(equivalent_labels_, labels_[i][j]);
      }
    }
  }

  // Перемаппинг меток
  std::vector<int> remap(label_count, 0);
  int current_label = 1;
  for (int i = 1; i < label_count; ++i) {
    if (equivalent_labels_[i] == i) {
      remap[i] = current_label;
      current_label++;
    }
  }

  // Применение перемаппинга
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (labels_[i][j] != 0) {
        labels_[i][j] = remap[labels_[i][j]];
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
      output.push_back(labels_[i][j]);
    }
  }

  return true;
}

}  // namespace artyushkina_markirovka
