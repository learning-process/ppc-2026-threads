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

std::vector<NeighborOffsetAll> GetFirstPassNeighbors() {
  std::vector<NeighborOffsetAll> neighbors(4);
  neighbors[0] = {.di = -1, .dj = -1, .check_i_min = 1, .check_i_max = 0, .check_j_min = 1, .check_j_max = 0};
  neighbors[1] = {.di = -1, .dj = 0,  .check_i_min = 1, .check_i_max = 0, .check_j_min = 0, .check_j_max = 0};
  neighbors[2] = {.di = -1, .dj = 1,  .check_i_min = 1, .check_i_max = 0, .check_j_min = 0, .check_j_max = 1};
  neighbors[3] = {.di = 0,  .dj = -1, .check_i_min = 0, .check_i_max = 0, .check_j_min = 1, .check_j_max = 0};
  return neighbors;
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

  if (labels_[static_cast<size_t>(ni)][static_cast<size_t>(nj)] != 0) {
    int label = labels_[static_cast<size_t>(ni)][static_cast<size_t>(nj)];
    neighbor_labels.push_back(label);
    min_label = std::min(label, min_label);
  }
}

void MarkingComponentsALL::FirstPass() {
  const auto &input = GetInput();
  int next_label = 1;
  std::vector<NeighborOffsetAll> first_pass_neighbors = GetFirstPassNeighbors();

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      size_t idx = (static_cast<size_t>(i) * static_cast<size_t>(cols_)) + static_cast<size_t>(j) + 2;

      if (input[idx] != 0) {
        continue;
      }
      
      std::vector<int> neighbor_labels;
      int min_label = next_label;

      for (const auto &neighbor : first_pass_neighbors) {
        ProcessNeighborFirstPass(i, j, neighbor, neighbor_labels, min_label);
      }

      if (neighbor_labels.empty()) {
        labels_[static_cast<size_t>(i)][static_cast<size_t>(j)] = next_label;
        equivalent_labels_.push_back(next_label);
        ++next_label;
      } else {
        labels_[static_cast<size_t>(i)][static_cast<size_t>(j)] = min_label;
        for (int label : neighbor_labels) {
          if (label != min_label) {
            UnionLabels(equivalent_labels_, label, min_label);
          }
        }
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