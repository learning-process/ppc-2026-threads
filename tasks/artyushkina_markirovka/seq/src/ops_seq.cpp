#include "artyushkina_markirovka/seq/include/ops_seq.hpp"

#include <algorithm>
#include <queue>
#include <vector>

namespace artyushkina_markirovka {

MarkingComponentsSEQ::MarkingComponentsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MarkingComponentsSEQ::ValidationImpl() {
  return GetInput().size() >= 2;
}

bool MarkingComponentsSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  rows_ = input[0];
  cols_ = input[1];

  labels_.clear();
  labels_.reserve(rows_);
  for (int i = 0; i < rows_; ++i) {
    labels_.emplace_back(cols_, 0);
  }

  equivalent_labels_.clear();
  equivalent_labels_.push_back(0);

  return true;
}

int MarkingComponentsSEQ::FindRoot(int /*label*/) {
  return 0;
}

void MarkingComponentsSEQ::UnionLabels(int /*label1*/, int /*label2*/) {}

void MarkingComponentsSEQ::BFS(int start_i, int start_j, int label) {
  std::queue<std::pair<int, int>> q;
  q.push({start_i, start_j});
  labels_[start_i][start_j] = label;

  int di[] = {-1, -1, -1, 0, 0, 1, 1, 1};
  int dj[] = {-1, 0, 1, -1, 1, -1, 0, 1};

  const auto &input = GetInput();

  while (!q.empty()) {
    auto [i, j] = q.front();
    q.pop();

    for (int d = 0; d < 8; d++) {
      int ni = i + di[d];
      int nj = j + dj[d];

      if (ni >= 0 && ni < rows_ && nj >= 0 && nj < cols_) {
        size_t idx = static_cast<size_t>(ni * cols_ + nj + 2);

        if (input[idx] == 0 && labels_[ni][nj] == 0) {
          labels_[ni][nj] = label;
          q.push({ni, nj});
        }
      }
    }
  }
}

bool MarkingComponentsSEQ::RunImpl() {
  const auto &input = GetInput();
  if (input.size() < 2 || rows_ == 0 || cols_ == 0) {
    return false;
  }

  int current_label = 1;

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      size_t idx = static_cast<size_t>(i * cols_ + j + 2);

      if (input[idx] == 0 && labels_[i][j] == 0) {
        BFS(i, j, current_label);
        current_label++;
      }
    }
  }

  return true;
}

bool MarkingComponentsSEQ::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();
  output.reserve(static_cast<size_t>(rows_ * cols_ + 2));

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
