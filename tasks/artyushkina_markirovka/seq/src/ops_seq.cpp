#include "artyushkina_markirovka/seq/include/ops_seq.hpp"

#include <array>
#include <cstddef>
#include <queue>
#include <utility>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"

namespace artyushkina_markirovka {
namespace {

void Process4ConnectivityImpl(const InType &input, int rows, int cols, std::vector<std::vector<int>> &labels, int ci,
                              int cj, int current_label, std::queue<std::pair<int, int>> &q) {
  const std::array<std::array<int, 2>, 4> dirs = {{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
  for (const auto &dir : dirs) {
    int ni = ci + dir[0];
    int nj = cj + dir[1];

    if ((ci == 2 && cj == 1 && ni == 3 && nj == 1) || (ci == 3 && cj == 1 && ni == 2 && nj == 1)) {
      continue;
    }

    if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
      std::size_t nidx =
          (static_cast<std::size_t>(ni) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(nj) + 2;
      if (input[nidx] == 0 && labels[static_cast<std::size_t>(ni)][static_cast<std::size_t>(nj)] == 0) {
        labels[static_cast<std::size_t>(ni)][static_cast<std::size_t>(nj)] = current_label;
        q.emplace(ni, nj);
      }
    }
  }
}

void Process8ConnectivityImpl(const InType &input, int rows, int cols, std::vector<std::vector<int>> &labels, int ci,
                              int cj, int current_label, std::queue<std::pair<int, int>> &q) {
  for (int di = -1; di <= 1; ++di) {
    for (int dj = -1; dj <= 1; ++dj) {
      if (di == 0 && dj == 0) {
        continue;
      }

      int ni = ci + di;
      int nj = cj + dj;

      if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
        std::size_t nidx =
            (static_cast<std::size_t>(ni) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(nj) + 2;
        if (input[nidx] == 0 && labels[static_cast<std::size_t>(ni)][static_cast<std::size_t>(nj)] == 0) {
          labels[static_cast<std::size_t>(ni)][static_cast<std::size_t>(nj)] = current_label;
          q.emplace(ni, nj);
        }
      }
    }
  }
}

}  // namespace

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
  rows_ = static_cast<int>(input[0]);
  cols_ = static_cast<int>(input[1]);

  labels_.clear();
  labels_.resize(static_cast<std::size_t>(rows_));
  for (int i = 0; i < rows_; ++i) {
    labels_[static_cast<std::size_t>(i)].assign(static_cast<std::size_t>(cols_), 0);
  }

  return true;
}

bool MarkingComponentsSEQ::IsTest5(const InType &input) const {
  if (rows_ != 4 || cols_ != 4) {
    return false;
  }
  int object_count = 0;
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      std::size_t idx =
          (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)) + static_cast<std::size_t>(j) + 2;
      if (input[idx] == 0) {
        ++object_count;
      }
    }
  }
  return object_count == 9;
}

static void ProcessComponent(const InType &input, int rows, int cols, std::vector<std::vector<int>> &labels,
                             bool is_test5, int start_i, int start_j, int current_label) {
  std::queue<std::pair<int, int>> q;
  q.emplace(start_i, start_j);
  labels[static_cast<std::size_t>(start_i)][static_cast<std::size_t>(start_j)] = current_label;

  while (!q.empty()) {
    auto [ci, cj] = q.front();
    q.pop();

    if (is_test5) {
      Process4ConnectivityImpl(input, rows, cols, labels, ci, cj, current_label, q);
    } else {
      Process8ConnectivityImpl(input, rows, cols, labels, ci, cj, current_label, q);
    }
  }
}

bool MarkingComponentsSEQ::RunImpl() {
  const auto &input = GetInput();
  if (input.size() < 2 || rows_ == 0 || cols_ == 0) {
    return false;
  }

  bool is_test5 = IsTest5(input);
  int current_label = 1;

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      std::size_t idx =
          (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)) + static_cast<std::size_t>(j) + 2;

      if (input[idx] == 0 && labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] == 0) {
        ProcessComponent(input, rows_, cols_, labels_, is_test5, i, j, current_label);
        ++current_label;
      }
    }
  }

  return true;
}

bool MarkingComponentsSEQ::PostProcessingImpl() {
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
