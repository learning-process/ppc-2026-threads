#include "gaivoronskiy_m_marking_binary_components/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <queue>
#include <utility>
#include <vector>

#include "gaivoronskiy_m_marking_binary_components/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gaivoronskiy_m_marking_binary_components {

GaivoronskiyMMarkingBinaryComponentsOMP::GaivoronskiyMMarkingBinaryComponentsOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GaivoronskiyMMarkingBinaryComponentsOMP::ValidationImpl() {
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

bool GaivoronskiyMMarkingBinaryComponentsOMP::PreProcessingImpl() {
  const auto &input = GetInput();
  int rows = input[0];
  int cols = input[1];
  GetOutput().assign((static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols)) + 2, 0);
  GetOutput()[0] = rows;
  GetOutput()[1] = cols;
  return true;
}

namespace {

constexpr std::array<int, 4> kDx = {-1, 1, 0, 0};
constexpr std::array<int, 4> kDy = {0, 0, -1, 1};

void BfsLabelInStrip(const InType &input, std::vector<int> &local_plane, int cols, int r_begin, int r_end,
                     int start_row, int start_col, int label) {
  std::queue<std::pair<int, int>> queue;
  queue.emplace(start_row, start_col);
  local_plane[(start_row * cols) + start_col] = label;

  while (!queue.empty()) {
    auto [cx, cy] = queue.front();
    queue.pop();

    for (std::size_t dir = 0; dir < 4; dir++) {
      int nx = cx + kDx.at(dir);
      int ny = cy + kDy.at(dir);
      if (nx < r_begin || nx >= r_end || ny < 0 || ny >= cols) {
        continue;
      }
      int n_flat = (nx * cols) + ny;
      int nidx = n_flat + 2;
      if (input[nidx] == 0 && local_plane[n_flat] == 0) {
        local_plane[n_flat] = label;
        queue.emplace(nx, ny);
      }
    }
  }
}

int FindRoot(std::vector<int> &parent, int x) {
  int root = x;
  while (parent[static_cast<std::size_t>(root)] != root) {
    root = parent[static_cast<std::size_t>(root)];
  }
  while (parent[static_cast<std::size_t>(x)] != x) {
    int p = parent[static_cast<std::size_t>(x)];
    parent[static_cast<std::size_t>(x)] = root;
    x = p;
  }
  return root;
}

void UniteLabels(std::vector<int> &parent, int a, int b) {
  a = FindRoot(parent, a);
  b = FindRoot(parent, b);
  if (a == b) {
    return;
  }
  if (a < b) {
    parent[static_cast<std::size_t>(b)] = a;
  } else {
    parent[static_cast<std::size_t>(a)] = b;
  }
}

}  // namespace

bool GaivoronskiyMMarkingBinaryComponentsOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  const int rows = input[0];
  const int cols = input[1];
  const int cells = rows * cols;

  if (cells == 0) {
    return true;
  }

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads < 1) {
    num_threads = 1;
  }
  num_threads = std::min(num_threads, rows);

  std::vector<int> row_starts(static_cast<std::size_t>(num_threads) + 1);
  for (int t = 0; t <= num_threads; t++) {
    row_starts[static_cast<std::size_t>(t)] = (t * rows) / num_threads;
  }

  std::vector<std::vector<int>> local_planes(static_cast<std::size_t>(num_threads),
                                             std::vector<int>(static_cast<std::size_t>(cells), 0));
  std::vector<int> labels_used(static_cast<std::size_t>(num_threads), 0);

#pragma omp parallel num_threads(num_threads) default(none) shared(input, local_planes, labels_used, row_starts, rows, \
                                                                     cols, num_threads)
  {
    const int tid = omp_get_thread_num();
    const int r_begin = row_starts[static_cast<std::size_t>(tid)];
    const int r_end = row_starts[static_cast<std::size_t>(tid) + 1];
    int next_label = 0;

    if (r_begin < r_end) {
      auto &plane = local_planes[static_cast<std::size_t>(tid)];
      for (int r = r_begin; r < r_end; r++) {
        for (int c = 0; c < cols; c++) {
          const int flat = (r * cols) + c;
          const int idx = flat + 2;
          if (input[idx] == 0 && plane[static_cast<std::size_t>(flat)] == 0) {
            next_label++;
            BfsLabelInStrip(input, plane, cols, r_begin, r_end, r, c, next_label);
          }
        }
      }
    }

    labels_used[static_cast<std::size_t>(tid)] = next_label;
  }

  std::vector<int> base(static_cast<std::size_t>(num_threads), 0);
  int sum = 0;
  for (int t = 0; t < num_threads; t++) {
    base[static_cast<std::size_t>(t)] = sum;
    sum += labels_used[static_cast<std::size_t>(t)];
  }

  const int max_global_before_merge = sum;

  if (max_global_before_merge == 0) {
    return true;
  }

  for (int t = 0; t < num_threads; t++) {
    const int r_begin = row_starts[static_cast<std::size_t>(t)];
    const int r_end = row_starts[static_cast<std::size_t>(t) + 1];
    if (r_begin >= r_end) {
      continue;
    }
    const auto &plane = local_planes[static_cast<std::size_t>(t)];
    const int offset = base[static_cast<std::size_t>(t)];
    for (int r = r_begin; r < r_end; r++) {
      for (int c = 0; c < cols; c++) {
        const int flat = (r * cols) + c;
        const int loc = plane[static_cast<std::size_t>(flat)];
        if (loc > 0) {
          output[flat + 2] = offset + loc;
        }
      }
    }
  }

  std::vector<int> parent(static_cast<std::size_t>(max_global_before_merge) + 1);
  for (int i = 1; i <= max_global_before_merge; i++) {
    parent[static_cast<std::size_t>(i)] = i;
  }

  for (int t = 0; t + 1 < num_threads; t++) {
    const int boundary_row = row_starts[static_cast<std::size_t>(t + 1)];
    if (boundary_row <= 0 || boundary_row >= rows) {
      continue;
    }
    const int ra = boundary_row - 1;
    const int rb = boundary_row;
    for (int c = 0; c < cols; c++) {
      const int ia = (ra * cols) + c + 2;
      const int ib = (rb * cols) + c + 2;
      if (input[ia] == 0 && input[ib] == 0) {
        const int la = output[ia];
        const int lb = output[ib];
        if (la > 0 && lb > 0) {
          UniteLabels(parent, la, lb);
        }
      }
    }
  }

  for (int i = 1; i <= max_global_before_merge; i++) {
    parent[static_cast<std::size_t>(i)] = FindRoot(parent, i);
  }

  std::vector<int> remap(static_cast<std::size_t>(max_global_before_merge) + 1, 0);
  int next_final = 1;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      const int idx = (r * cols) + c + 2;
      if (input[idx] != 0) {
        continue;
      }
      const int lbl = output[idx];
      if (lbl <= 0) {
        continue;
      }
      const int root = parent[static_cast<std::size_t>(lbl)];
      if (remap[static_cast<std::size_t>(root)] == 0) {
        remap[static_cast<std::size_t>(root)] = next_final++;
      }
    }
  }

#pragma omp parallel for default(none) shared(input, output, parent, remap, rows, cols) schedule(static)
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      const int idx = (r * cols) + c + 2;
      if (input[idx] != 0) {
        output[idx] = 0;
      } else {
        const int lbl = output[idx];
        const int root = parent[static_cast<std::size_t>(lbl)];
        output[idx] = remap[static_cast<std::size_t>(root)];
      }
    }
  }

  return true;
}

bool GaivoronskiyMMarkingBinaryComponentsOMP::PostProcessingImpl() {
  return true;
}

}  // namespace gaivoronskiy_m_marking_binary_components
