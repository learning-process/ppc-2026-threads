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

int NumThreadsForRows(int rows) {
  const int capped = std::max(ppc::util::GetNumThreads(), 1);
  return std::min(capped, rows);
}

std::vector<int> MakeRowStarts(int rows, int num_threads) {
  std::vector<int> row_starts(static_cast<std::size_t>(num_threads) + 1);
  for (int thread_idx = 0; thread_idx <= num_threads; thread_idx++) {
    row_starts[static_cast<std::size_t>(thread_idx)] = (thread_idx * rows) / num_threads;
  }
  return row_starts;
}

void ParallelLabelStrips(const InType &input, int cols, int num_threads, const std::vector<int> &row_starts,
                         std::vector<std::vector<int>> &local_planes, std::vector<int> &labels_used) {
#pragma omp parallel num_threads(num_threads) default(none) \
    shared(input, local_planes, labels_used, row_starts, cols, num_threads)
  {
    const int tid = omp_get_thread_num();
    const int r_begin = row_starts[static_cast<std::size_t>(tid)];
    const int r_end = row_starts[static_cast<std::size_t>(tid) + 1];
    int next_label = 0;

    if (r_begin < r_end) {
      auto &plane = local_planes[static_cast<std::size_t>(tid)];
      for (int row = r_begin; row < r_end; row++) {
        for (int col = 0; col < cols; col++) {
          const int flat = (row * cols) + col;
          const int idx = flat + 2;
          if (input[idx] == 0 && plane[static_cast<std::size_t>(flat)] == 0) {
            next_label++;
            BfsLabelInStrip(input, plane, cols, r_begin, r_end, row, col, next_label);
          }
        }
      }
    }

    labels_used[static_cast<std::size_t>(tid)] = next_label;
  }
}

std::pair<std::vector<int>, int> BuildLabelBases(const std::vector<int> &labels_used, int num_threads) {
  std::vector<int> base(static_cast<std::size_t>(num_threads), 0);
  int sum = 0;
  for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
    base[static_cast<std::size_t>(thread_idx)] = sum;
    sum += labels_used[static_cast<std::size_t>(thread_idx)];
  }
  return {base, sum};
}

void CopyStripsToGlobalOutput(OutType &output, int cols, int num_threads, const std::vector<int> &row_starts,
                              const std::vector<int> &base, const std::vector<std::vector<int>> &local_planes) {
  for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
    const int r_begin = row_starts[static_cast<std::size_t>(thread_idx)];
    const int r_end = row_starts[static_cast<std::size_t>(thread_idx) + 1U];
    if (r_begin >= r_end) {
      continue;
    }
    const auto &plane = local_planes[static_cast<std::size_t>(thread_idx)];
    const int offset = base[static_cast<std::size_t>(thread_idx)];
    for (int row = r_begin; row < r_end; row++) {
      for (int col = 0; col < cols; col++) {
        const int flat = (row * cols) + col;
        const int loc = plane[static_cast<std::size_t>(flat)];
        if (loc > 0) {
          output[flat + 2] = offset + loc;
        }
      }
    }
  }
}

void MergeBoundariesUnionFind(const InType &input, OutType &output, int rows, int cols, int num_threads,
                              const std::vector<int> &row_starts, std::vector<int> &parent) {
  for (int thread_idx = 0; thread_idx + 1 < num_threads; thread_idx++) {
    const int next_strip = thread_idx + 1;
    const int boundary_row = row_starts[static_cast<std::size_t>(next_strip)];
    if (boundary_row <= 0 || boundary_row >= rows) {
      continue;
    }
    const int ra = boundary_row - 1;
    const int rb = boundary_row;
    for (int col = 0; col < cols; col++) {
      const int ia = (ra * cols) + col + 2;
      const int ib = (rb * cols) + col + 2;
      if (input[ia] != 0 || input[ib] != 0) {
        continue;
      }
      const int la = output[ia];
      const int lb = output[ib];
      if (la > 0 && lb > 0) {
        UniteLabels(parent, la, lb);
      }
    }
  }
}

void FlattenParentForest(std::vector<int> &parent, int max_label) {
  for (int label_idx = 1; label_idx <= max_label; label_idx++) {
    parent[static_cast<std::size_t>(label_idx)] = FindRoot(parent, label_idx);
  }
}

std::vector<int> BuildFinalRemap(const InType &input, const OutType &output, const std::vector<int> &parent, int rows,
                                 int cols, int max_label) {
  std::vector<int> remap(static_cast<std::size_t>(max_label) + 1, 0);
  int next_final = 1;
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      const int idx = (row * cols) + col + 2;
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
  return remap;
}

void ApplyRemapInParallel(const InType &input, OutType &output, const std::vector<int> &parent,
                          const std::vector<int> &remap, int rows, int cols) {
#pragma omp parallel for default(none) shared(input, output, parent, remap, rows, cols) schedule(static)
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      const int idx = (row * cols) + col + 2;
      if (input[idx] != 0) {
        output[idx] = 0;
      } else {
        const int lbl = output[idx];
        const int root = parent[static_cast<std::size_t>(lbl)];
        output[idx] = remap[static_cast<std::size_t>(root)];
      }
    }
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

  const int num_threads = NumThreadsForRows(rows);
  const std::vector<int> row_starts = MakeRowStarts(rows, num_threads);

  std::vector<std::vector<int>> local_planes(static_cast<std::size_t>(num_threads),
                                             std::vector<int>(static_cast<std::size_t>(cells), 0));
  std::vector<int> labels_used(static_cast<std::size_t>(num_threads), 0);

  ParallelLabelStrips(input, cols, num_threads, row_starts, local_planes, labels_used);

  const auto [base, max_global_before_merge] = BuildLabelBases(labels_used, num_threads);
  if (max_global_before_merge == 0) {
    return true;
  }

  CopyStripsToGlobalOutput(output, cols, num_threads, row_starts, base, local_planes);

  std::vector<int> parent(static_cast<std::size_t>(max_global_before_merge) + 1);
  for (int label_idx = 1; label_idx <= max_global_before_merge; label_idx++) {
    parent[static_cast<std::size_t>(label_idx)] = label_idx;
  }

  MergeBoundariesUnionFind(input, output, rows, cols, num_threads, row_starts, parent);
  FlattenParentForest(parent, max_global_before_merge);

  const std::vector<int> remap = BuildFinalRemap(input, output, parent, rows, cols, max_global_before_merge);
  ApplyRemapInParallel(input, output, parent, remap, rows, cols);

  return true;
}

bool GaivoronskiyMMarkingBinaryComponentsOMP::PostProcessingImpl() {
  return true;
}

}  // namespace gaivoronskiy_m_marking_binary_components
