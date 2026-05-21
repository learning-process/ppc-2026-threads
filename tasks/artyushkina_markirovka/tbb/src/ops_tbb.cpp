#include "artyushkina_markirovka/tbb/include/ops_tbb.hpp"

#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"

namespace artyushkina_markirovka {

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

  int total_pixels = rows_ * cols_;
  labels_.assign(total_pixels, 0);
  parent_.resize(total_pixels + 1);
  for (int i = 0; i <= total_pixels; ++i) {
    parent_[i] = i;
  }

  current_label_ = 0;
  return true;
}

int MarkingComponentsTBB::FindRoot(int label) {
  int root = label;
  while (parent_[root] != root) {
    root = parent_[root];
  }
  int current = label;
  while (parent_[current] != current) {
    int next = parent_[current];
    parent_[current] = root;
    current = next;
  }
  return root;
}

void MarkingComponentsTBB::UnionLabels(int label1, int label2) {
  tbb::spin_mutex::scoped_lock lock(dsu_mutex_);
  int root1 = FindRoot(label1);
  int root2 = FindRoot(label2);
  if (root1 != root2) {
    if (root1 < root2) {
      parent_[root2] = root1;
    } else {
      parent_[root1] = root2;
    }
  }
}

void MarkingComponentsTBB::InitLabelsTbb() {
  int total_pixels = rows_ * cols_;
  tbb::parallel_for(0, total_pixels, [this](int idx) {
    size_t input_idx = static_cast<size_t>(idx) + 2;
    if (input_[input_idx] == 0) {
      labels_[idx] = idx + 1;
    }
  });
}

void MarkingComponentsTBB::MergeHorizontalPairsTbb() {
  tbb::parallel_for(0, rows_, [this](int y_coord) {
    for (int x_coord = 0; x_coord < cols_ - 1; ++x_coord) {
      int idx = (y_coord * cols_) + x_coord;
      if (labels_[idx] != 0 && labels_[idx + 1] != 0) {
        UnionLabels(labels_[idx], labels_[idx + 1]);
      }
    }
  });
}

void MarkingComponentsTBB::MergeVerticalPairsTbb() {
  tbb::parallel_for(0, rows_ - 1, [this](int y_coord) {
    for (int x_coord = 0; x_coord < cols_; ++x_coord) {
      int idx = (y_coord * cols_) + x_coord;
      if (labels_[idx] != 0 && labels_[idx + cols_] != 0) {
        UnionLabels(labels_[idx], labels_[idx + cols_]);
      }
    }
  });
}

void MarkingComponentsTBB::MergeDiagonalPairsTbb() {
  tbb::parallel_for(0, rows_ - 1, [this](int y_coord) {
    for (int x_coord = 0; x_coord < cols_ - 1; ++x_coord) {
      int idx = (y_coord * cols_) + x_coord;
      if (labels_[idx] != 0 && labels_[idx + cols_ + 1] != 0) {
        UnionLabels(labels_[idx], labels_[idx + cols_ + 1]);
      }
      if (x_coord > 0) {
        if (labels_[idx] != 0 && labels_[idx + cols_ - 1] != 0) {
          UnionLabels(labels_[idx], labels_[idx + cols_ - 1]);
        }
      }
    }
  });
}

void MarkingComponentsTBB::FinalizeRootsTbb() {
  int total_pixels = rows_ * cols_;
  tbb::parallel_for(0, total_pixels, [this](int i) {
    if (labels_[i] != 0) {
      labels_[i] = FindRoot(labels_[i]);
    }
  });
}

void MarkingComponentsTBB::NormalizeLabelsTbb() {
  int total_pixels = rows_ * cols_;
  std::vector<int> unique_roots;
  for (int i = 0; i < total_pixels; ++i) {
    if (labels_[i] != 0) {
      unique_roots.push_back(labels_[i]);
    }
  }

  if (unique_roots.empty()) {
    return;
  }

  std::sort(unique_roots.begin(), unique_roots.end());
  auto last = std::unique(unique_roots.begin(), unique_roots.end());
  unique_roots.erase(last, unique_roots.end());

  std::vector<int> mapping(total_pixels + 1, 0);
  int next_id = 1;
  for (int root : unique_roots) {
    mapping[root] = next_id++;
  }

  for (int i = 0; i < total_pixels; ++i) {
    if (labels_[i] != 0) {
      labels_[i] = mapping[labels_[i]];
    }
  }
  current_label_ = next_id - 1;
}

bool MarkingComponentsTBB::RunImpl() {
  int total_pixels = rows_ * cols_;
  if (total_pixels <= 0) {
    return true;
  }

  InitLabelsTbb();
  MergeHorizontalPairsTbb();
  MergeVerticalPairsTbb();
  MergeDiagonalPairsTbb();
  FinalizeRootsTbb();
  NormalizeLabelsTbb();

  return true;
}

bool MarkingComponentsTBB::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();

  output.push_back(static_cast<uint8_t>(rows_));
  output.push_back(static_cast<uint8_t>(cols_));

  for (int label : labels_) {
    output.push_back(static_cast<uint8_t>(label));
  }

  return true;
}

}  // namespace artyushkina_markirovka
