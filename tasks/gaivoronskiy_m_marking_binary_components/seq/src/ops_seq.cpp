#include "gaivoronskiy_m_marking_binary_components/seq/include/ops_seq.hpp"

#include <queue>
#include <utility>
#include <vector>

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
  return static_cast<int>(input.size()) == rows * cols + 2;
}

bool GaivoronskiyMMarkingBinaryComponentsSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  int rows = input[0];
  int cols = input[1];
  GetOutput().assign(static_cast<size_t>(rows * cols + 2), 0);
  GetOutput()[0] = rows;
  GetOutput()[1] = cols;
  return true;
}

bool GaivoronskiyMMarkingBinaryComponentsSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  int rows = input[0];
  int cols = input[1];

  static constexpr int kDx[] = {-1, 1, 0, 0};
  static constexpr int kDy[] = {0, 0, -1, 1};

  int label = 0;

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      int idx = i * cols + j + 2;
      if (input[idx] == 0 && output[idx] == 0) {
        label++;
        std::queue<std::pair<int, int>> q;
        q.emplace(i, j);
        output[idx] = label;

        while (!q.empty()) {
          auto [x, y] = q.front();
          q.pop();

          for (int d = 0; d < 4; d++) {
            int nx = x + kDx[d];
            int ny = y + kDy[d];
            if (nx >= 0 && nx < rows && ny >= 0 && ny < cols) {
              int nidx = nx * cols + ny + 2;
              if (input[nidx] == 0 && output[nidx] == 0) {
                output[nidx] = label;
                q.emplace(nx, ny);
              }
            }
          }
        }
      }
    }
  }

  return true;
}

bool GaivoronskiyMMarkingBinaryComponentsSEQ::PostProcessingImpl() { return true; }

}  // namespace gaivoronskiy_m_marking_binary_components
