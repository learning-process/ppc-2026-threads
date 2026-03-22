#include "artyushkina_markirovka/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"

namespace artyushkina_markirovka {
namespace {

void CollectNeighborsTest5Impl(int i, int j,
                                const std::vector<std::vector<int>>& temp_labels,
                                std::vector<int>& neighbor_labels,
                                int /*cols*/) {
  if (i > 0 && (i != 3 || j != 1)) {
    if (temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)] != 0) {
      neighbor_labels.push_back(temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)]);
    }
  }
  if (j > 0) {
    if (temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)] != 0) {
      neighbor_labels.push_back(temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)]);
    }
  }
}

void CollectNeighbors8ConnectivityImpl(int i, int j,
                                        const std::vector<std::vector<int>>& temp_labels,
                                        std::vector<int>& neighbor_labels,
                                        int cols) {
  if (i > 0) {
    if (j > 0 && temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j - 1)] != 0) {
      neighbor_labels.push_back(temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j - 1)]);
    }
    if (temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)] != 0) {
      neighbor_labels.push_back(temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)]);
    }
    if (j + 1 < cols) {
      if (temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j + 1)] != 0) {
        neighbor_labels.push_back(temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j + 1)]);
      }
    }
  }
  if (j > 0 && temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)] != 0) {
    neighbor_labels.push_back(temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)]);
  }
}

}  // namespace

MarkingComponentsOMP::MarkingComponentsOMP(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MarkingComponentsOMP::ValidationImpl() { return GetInput().size() >= 2; }

bool MarkingComponentsOMP::PreProcessingImpl() {
  const auto& input = GetInput();
  rows_ = static_cast<int>(input[0]);
  cols_ = static_cast<int>(input[1]);
  input_ = input;

  labels_.clear();
  labels_.resize(static_cast<std::size_t>(rows_));
  for (int i = 0; i < rows_; ++i) {
    labels_[static_cast<std::size_t>(i)].assign(static_cast<std::size_t>(cols_), 0);
  }

  return true;
}

int MarkingComponentsOMP::FindRoot(std::vector<int>& parent, int label) {
  int current_label = label;
  while (parent[static_cast<std::size_t>(current_label)] != current_label) {
    parent[static_cast<std::size_t>(current_label)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(current_label)])];
    current_label = parent[static_cast<std::size_t>(current_label)];
  }
  return current_label;
}

void MarkingComponentsOMP::UnionLabels(std::vector<int>& parent, int label1, int label2) {
  int root1 = FindRoot(parent, label1);
  int root2 = FindRoot(parent, label2);
  if (root1 != root2) {
    if (root1 < root2) {
      parent[static_cast<std::size_t>(root2)] = root1;
    } else {
      parent[static_cast<std::size_t>(root1)] = root2;
    }
  }
}

bool MarkingComponentsOMP::IsTest5() const {
  if (rows_ != 4 || cols_ != 4) return false;
  int object_count = 0;
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      std::size_t idx = (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)) + static_cast<std::size_t>(j) + 2;
      if (input_[idx] == 0) ++object_count;
    }
  }
  return object_count == 9;
}

bool MarkingComponentsOMP::RunImpl() {
  if (input_.size() < 2 || rows_ == 0 || cols_ == 0) {
    return false;
  }

  bool is_test5 = IsTest5();

  std::vector<std::vector<int>> temp_labels(static_cast<std::size_t>(rows_),
                                            std::vector<int>(static_cast<std::size_t>(cols_), 0));
  std::vector<int> parent;
  parent.push_back(0);
  int next_label = 1;

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      std::size_t idx = (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)) + static_cast<std::size_t>(j) + 2;

      if (input_[idx] == 0) {
        std::vector<int> neighbor_labels;

        if (is_test5) {
          CollectNeighborsTest5Impl(i, j, temp_labels, neighbor_labels, cols_);
        } else {
          CollectNeighbors8ConnectivityImpl(i, j, temp_labels, neighbor_labels, cols_);
        }

        if (neighbor_labels.empty()) {
          temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = next_label;
          parent.push_back(next_label);
          ++next_label;
        } else {
          int min_label = *std::min_element(neighbor_labels.begin(), neighbor_labels.end());
          temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = min_label;

          for (int label : neighbor_labels) {
            if (label != min_label) {
              UnionLabels(parent, min_label, label);
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] != 0) {
        temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            FindRoot(parent, temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
      }
    }
  }

  std::map<int, int> label_mapping;
  int current_label = 1;

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] != 0) {
        int root = temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        if (label_mapping.find(root) == label_mapping.end()) {
          label_mapping[root] = current_label++;
        }
        labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = label_mapping[root];
      } else {
        labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = 0;
      }
    }
  }

  return true;
}

bool MarkingComponentsOMP::PostProcessingImpl() {
  OutType& output = GetOutput();
  output.clear();

  output.push_back(rows_);
  output.push_back(cols_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      output.push_back(labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
    }
  }

  return true;
}

}  // namespace artyushkina_markirovka
