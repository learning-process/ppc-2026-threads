#include "artyushkina_markirovka/seq/include/ops_seq.hpp"

#include <array>
#include <cstddef>
#include <queue>
#include <utility>

#include "artyushkina_markirovka/common/include/common.hpp"

namespace artyushkina_markirovka {

struct NeighborOffset {
  int di;
  int dj;
  bool check_i_min;
  bool check_i_max;
  bool check_j_min;
  bool check_j_max;
};

const std::array<NeighborOffset, 8> kNeighbors = {
    {NeighborOffset{
         .di = -1, .dj = -1, .check_i_min = true, .check_i_max = false, .check_j_min = true, .check_j_max = false},
     NeighborOffset{
         .di = -1, .dj = 0, .check_i_min = true, .check_i_max = false, .check_j_min = false, .check_j_max = false},
     NeighborOffset{
         .di = -1, .dj = 1, .check_i_min = true, .check_i_max = false, .check_j_min = false, .check_j_max = true},
     NeighborOffset{
         .di = 0, .dj = -1, .check_i_min = false, .check_i_max = false, .check_j_min = true, .check_j_max = false},
     NeighborOffset{
         .di = 0, .dj = 1, .check_i_min = false, .check_i_max = false, .check_j_min = false, .check_j_max = true},
     NeighborOffset{
         .di = 1, .dj = -1, .check_i_min = false, .check_i_max = true, .check_j_min = true, .check_j_max = false},
     NeighborOffset{
         .di = 1, .dj = 0, .check_i_min = false, .check_i_max = true, .check_j_min = false, .check_j_max = false},
     NeighborOffset{
         .di = 1, .dj = 1, .check_i_min = false, .check_i_max = true, .check_j_min = false, .check_j_max = true}}};

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

bool MarkingComponentsSEQ::IsValidNeighbor(int i, int j, const NeighborOffset &offset) const {
  if (offset.check_i_min && i <= 0) {
    return false;
  }
  if (offset.check_i_max && i >= rows_ - 1) {
    return false;
  }
  if (offset.check_j_min && j <= 0) {
    return false;
  }
  if (offset.check_j_max && j >= cols_ - 1) {
    return false;
  }
  return true;
}

void MarkingComponentsSEQ::ProcessNeighbor(int i, int j, const NeighborOffset &offset, int label,
                                           std::queue<std::pair<int, int>> &q) {
  if (!IsValidNeighbor(i, j, offset)) {
    return;
  }

  int ni = i + offset.di;
  int nj = j + offset.dj;

  const auto &input = GetInput();
  std::size_t idx = (static_cast<std::size_t>(ni) * static_cast<std::size_t>(cols_)) + static_cast<std::size_t>(nj) + 2;

  if (input[idx] == 0 && labels_[ni][nj] == 0) {
    labels_[ni][nj] = label;
    q.emplace(ni, nj);
  }
}

void MarkingComponentsSEQ::BFS(int start_i, int start_j, int label) {
  std::queue<std::pair<int, int>> q;
  q.emplace(start_i, start_j);
  labels_[start_i][start_j] = label;

  while (!q.empty()) {
    auto [i, j] = q.front();
    q.pop();

    for (const auto &offset : kNeighbors) {
      ProcessNeighbor(i, j, offset, label, q);
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
      std::size_t idx =
          (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)) + static_cast<std::size_t>(j) + 2;

      if (input[idx] == 0 && labels_[i][j] == 0) {
        BFS(i, j, current_label);
        ++current_label;
      }
    }
  }

  return true;
}

bool MarkingComponentsSEQ::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();

  output.reserve((static_cast<std::size_t>(rows_) * static_cast<std::size_t>(cols_)) + 2);

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
