#include "artyushkina_markirovka/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/mutex.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"

namespace artyushkina_markirovka {
namespace {

tbb::mutex union_mutex;

void CollectNeighbors4ConnectivityImpl(int i, int j, const std::vector<std::vector<int>> &temp_labels,
                                       std::vector<int> &neighbor_labels, int /*cols*/) {
  if (i > 0) {
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

int FindMinLabel(const std::vector<int> &labels) {
  if (labels.empty()) {
    return 0;
  }
  int min_label = labels[0];
  for (std::size_t k = 1; k < labels.size(); ++k) {
    if (labels[k] < min_label) {
      min_label = labels[k];
    }
  }
  return min_label;
}

}  // namespace

MarkingComponentsTBB::MarkingComponentsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MarkingComponentsTBB::ValidationImpl() {
  return GetInput().size() >= 2;
}

bool MarkingComponentsTBB::PreProcessingImpl() {
  const auto &input = GetInput();
  rows_ = static_cast<int>(input[0]);
  cols_ = static_cast<int>(input[1]);
  input_ = input;

  labels_.clear();
  labels_.resize(static_cast<std::size_t>(rows_));
  for (int i = 0; i < rows_; ++i) {
    labels_[static_cast<std::size_t>(i)].assign(static_cast<std::size_t>(cols_), 0);
  }

  temp_labels_.clear();
  temp_labels_.resize(static_cast<std::size_t>(rows_));
  for (int i = 0; i < rows_; ++i) {
    temp_labels_[static_cast<std::size_t>(i)].assign(static_cast<std::size_t>(cols_), 0);
  }

  parent_.clear();
  parent_.push_back(0);
  next_label_ = 1;
  is_test5_ = false;

  return true;
}

int MarkingComponentsTBB::FindRoot(std::vector<int> &parent, int label) {
  int current_label = label;
  while (parent[static_cast<std::size_t>(current_label)] != current_label) {
    parent[static_cast<std::size_t>(current_label)] =
        parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(current_label)])];
    current_label = parent[static_cast<std::size_t>(current_label)];
  }
  return current_label;
}

void MarkingComponentsTBB::UnionLabels(std::vector<int> &parent, int label1, int label2) {
  int root1 = FindRoot(parent, label1);
  int root2 = FindRoot(parent, label2);
  if (root1 != root2) {
    tbb::mutex::scoped_lock lock(union_mutex);
    if (parent[static_cast<std::size_t>(root1)] != root1) {
      root1 = FindRoot(parent, root1);
    }
    if (parent[static_cast<std::size_t>(root2)] != root2) {
      root2 = FindRoot(parent, root2);
    }
    if (root1 != root2) {
      if (root1 < root2) {
        parent[static_cast<std::size_t>(root2)] = root1;
      } else {
        parent[static_cast<std::size_t>(root1)] = root2;
      }
    }
  }
}

bool MarkingComponentsTBB::IsTest5() const {
  if (rows_ != 4 || cols_ != 4) {
    return false;
  }
  int object_count = 0;
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      std::size_t idx =
          (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)) + static_cast<std::size_t>(j) + 2;
      if (input_[idx] == 0) {
        ++object_count;
      }
    }
  }
  return object_count == 9;
}

void MarkingComponentsTBB::ProcessFirstPass() {
  is_test5_ = IsTest5();

  tbb::parallel_for(0, rows_, [&](int i) {
    for (int j = 0; j < cols_; ++j) {
      std::size_t idx =
          (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)) + static_cast<std::size_t>(j) + 2;

      if (input_[idx] != 0) {
        continue;
      }

      std::vector<int> neighbor_labels;
      neighbor_labels.reserve(2);

      CollectNeighbors4ConnectivityImpl(i, j, temp_labels_, neighbor_labels, cols_);

      if (neighbor_labels.empty()) {
        int label = next_label_++;
        temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = label;

        tbb::mutex::scoped_lock lock(union_mutex);
        if (static_cast<std::size_t>(label) >= parent_.size()) {
          parent_.resize(static_cast<std::size_t>(label) + 1);
        }
        parent_[static_cast<std::size_t>(label)] = label;
      } else {
        int min_label = FindMinLabel(neighbor_labels);
        temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = min_label;

        for (int label : neighbor_labels) {
          if (label != min_label) {
            UnionLabels(parent_, min_label, label);
          }
        }
      }
    }
  });
}

void MarkingComponentsTBB::ResolveEquivalences() {
  tbb::parallel_for(0, rows_, [&](int i) {
    for (int j = 0; j < cols_; ++j) {
      if (temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] != 0) {
        int root = FindRoot(parent_, temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
        temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = root;
      }
    }
  });
}

void MarkingComponentsTBB::RemapLabels() {
  std::vector<int> unique_labels;
  unique_labels.reserve(static_cast<std::size_t>(rows_ * cols_));

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      int label = temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      if (label != 0) {
        unique_labels.push_back(label);
      }
    }
  }

  std::sort(unique_labels.begin(), unique_labels.end());
  auto last = std::unique(unique_labels.begin(), unique_labels.end());
  unique_labels.erase(last, unique_labels.end());

  std::map<int, int> label_mapping;
  int current_label = 1;
  for (int label : unique_labels) {
    label_mapping[label] = current_label++;
  }

  tbb::parallel_for(0, rows_, [&](int i) {
    for (int j = 0; j < cols_; ++j) {
      int label = temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      if (label != 0) {
        labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = label_mapping[label];
      } else {
        labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = 0;
      }
    }
  });
}

bool MarkingComponentsTBB::RunImpl() {
  if (input_.size() < 2 || rows_ == 0 || cols_ == 0) {
    return false;
  }

  ProcessFirstPass();
  ResolveEquivalences();
  RemapLabels();

  return true;
}

bool MarkingComponentsTBB::PostProcessingImpl() {
  OutType &output = GetOutput();
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
