#include "gaivoronskiy_m_marking_binary_components/seq/include/ops_seq.hpp"

#include <array>
#include <cstddef>
#include <queue>
#include <utility>
#include <vector>

#include "gaivoronskiy_m_marking_binary_components/common/include/common.hpp"

namespace gaivoronskiy_m_marking_binary_components {

GaivoronskiyMMarkingBinaryComponentsSEQ::GaivoronskiyMMarkingBinaryComponentsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GaivoronskiyMMarkingBinaryComponentsSEQ::ValidationImpl() {
  const auto &input = GetInput();
  if (input.size() < 2) {
    return false;
  }
  int rows = input[0];
  int cols = input[1];
  if (rows <= 0 || cols <= 0) {
    return false;
  }
  return static_cast<int>(input.size()) == (rows * cols) + 2;
}

bool GaivoronskiyMMarkingBinaryComponentsSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  int rows = input[0];
  int cols = input[1];
  GetOutput().assign((static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols)) + 2, 0);
  GetOutput()[0] = rows;
  GetOutput()[1] = cols;
  return true;
}

namespace {

void BfsLabel(const InType &input, OutType &output, int rows, int cols, int start_row, int start_col, int label) {
  static constexpr std::array<int, 4> kDx = {-1, 1, 0, 0};
  static constexpr std::array<int, 4> kDy = {0, 0, -1, 1};

  std::queue<std::pair<int, int>> queue;
  queue.emplace(start_row, start_col);
  output[(start_row * cols) + start_col + 2] = label;

  while (!queue.empty()) {
    auto [cx, cy] = queue.front();
    queue.pop();

    for (std::size_t dir = 0; dir < 4; dir++) {
      int nx = cx + kDx.at(dir);
      int ny = cy + kDy.at(dir);
      if (nx >= 0 && nx < rows && ny >= 0 && ny < cols) {
        int nidx = (nx * cols) + ny + 2;
        if (input[nidx] == 0 && output[nidx] == 0) {
          output[nidx] = label;
          queue.emplace(nx, ny);
        }
      }
    }
  }
}

}  // namespace

bool GaivoronskiyMMarkingBinaryComponentsSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  int rows = input[0];
  int cols = input[1];

  int label = 0;

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      int idx = (i * cols) + j + 2;
      if (input[idx] == 0 && output[idx] == 0) {
        label++;
        BfsLabel(input, output, rows, cols, i, j, label);
      }
    }
  }

  return true;
}

bool GaivoronskiyMMarkingBinaryComponentsSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace gaivoronskiy_m_marking_binary_components
