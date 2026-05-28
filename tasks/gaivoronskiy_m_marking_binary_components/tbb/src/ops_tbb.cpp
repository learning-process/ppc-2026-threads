#include "gaivoronskiy_m_marking_binary_components/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

#include "gaivoronskiy_m_marking_binary_components/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gaivoronskiy_m_marking_binary_components {

GaivoronskiyMMarkingBinaryComponentsTBB::GaivoronskiyMMarkingBinaryComponentsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GaivoronskiyMMarkingBinaryComponentsTBB::ValidationImpl() {
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

bool GaivoronskiyMMarkingBinaryComponentsTBB::PreProcessingImpl() {
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

void BfsLabelFull(const InType &input, OutType &output, int rows, int cols, int start_row, int start_col, int label) {
  std::vector<std::pair<int, int>> queue;
  queue.reserve(64);
  std::size_t head = 0;
  queue.emplace_back(start_row, start_col);
  output[(start_row * cols) + start_col + 2] = label;

  while (head < queue.size()) {
    const auto [cx, cy] = queue[head++];
    for (std::size_t dir = 0; dir < 4; dir++) {
      const int nx = cx + kDx.at(dir);
      const int ny = cy + kDy.at(dir);
      if (nx < 0 || nx >= rows || ny < 0 || ny >= cols) {
        continue;
      }
      const int nidx = (nx * cols) + ny + 2;
      if (input[nidx] == 0 && output[nidx] == 0) {
        output[nidx] = label;
        queue.emplace_back(nx, ny);
      }
    }
  }
}

void BfsLabelInStrip(const InType &input, OutType &output, int cols, int r_begin, int r_end, int start_row,
                     int start_col, int label) {
  std::vector<std::pair<int, int>> queue;
  queue.reserve(64);
  std::size_t head = 0;
  queue.emplace_back(start_row, start_col);
  output[(start_row * cols) + start_col + 2] = label;

  while (head < queue.size()) {
    const auto [cx, cy] = queue[head++];
    for (std::size_t dir = 0; dir < 4; dir++) {
      const int nx = cx + kDx.at(dir);
      const int ny = cy + kDy.at(dir);
      if (nx < r_begin || nx >= r_end || ny < 0 || ny >= cols) {
        continue;
      }
      const int nidx = (nx * cols) + ny + 2;
      if (input[nidx] == 0 && output[nidx] == 0) {
        output[nidx] = label;
        queue.emplace_back(nx, ny);
      }
    }
  }
}

void RunSequentialMarking(const InType &input, OutType &output, int rows, int cols) {
  int label = 0;
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      const int idx = (row * cols) + col + 2;
      if (input[idx] != 0 || output[idx] != 0) {
        continue;
      }
      label++;
      BfsLabelFull(input, output, rows, cols, row, col, label);
    }
  }
}

int FindRoot(std::vector<int> &parent, int x) {
  int root = x;
  while (parent[static_cast<std::size_t>(root)] != root) {
    root = parent[static_cast<std::size_t>(root)];
  }
  while (parent[static_cast<std::size_t>(x)] != x) {
    const int p = parent[static_cast<std::size_t>(x)];
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

int LabelSingleStrip(const InType &input, OutType &output, int cols, int r_begin, int r_end) {
  if (r_begin >= r_end) {
    return 0;
  }

  int next_label = 0;
  for (int row = r_begin; row < r_end; row++) {
    for (int col = 0; col < cols; col++) {
      const int idx = (row * cols) + col + 2;
      if (input[idx] != 0 || output[idx] != 0) {
        continue;
      }
      next_label++;
      BfsLabelInStrip(input, output, cols, r_begin, r_end, row, col, next_label);
    }
  }
  return next_label;
}

void LabelStripsRange(const InType &input, OutType &output, int cols, const std::vector<int> &row_starts,
                      std::vector<int> &labels_used, int begin_tid, int end_tid) {
  for (int tid = begin_tid; tid < end_tid; ++tid) {
    const int r_begin = row_starts[static_cast<std::size_t>(tid)];
    const int r_end = row_starts[static_cast<std::size_t>(tid) + 1];
    labels_used[static_cast<std::size_t>(tid)] = LabelSingleStrip(input, output, cols, r_begin, r_end);
  }
}

void ParallelLabelStrips(const InType &input, OutType &output, int cols, int num_threads,
                         const std::vector<int> &row_starts, std::vector<int> &labels_used) {
  tbb::task_arena arena(num_threads);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<int>(0, num_threads), [&](const tbb::blocked_range<int> &range) {
      LabelStripsRange(input, output, cols, row_starts, labels_used, range.begin(), range.end());
    });
  });
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

void AddStripBaseOffset(int thread_idx, OutType &output, int cols, const std::vector<int> &row_starts,
                        const std::vector<int> &base) {
  const int offset = base[static_cast<std::size_t>(thread_idx)];
  if (offset == 0) {
    return;
  }
  const int r_begin = row_starts[static_cast<std::size_t>(thread_idx)];
  const int r_end = row_starts[static_cast<std::size_t>(thread_idx) + 1U];
  for (int row = r_begin; row < r_end; row++) {
    for (int col = 0; col < cols; col++) {
      const int idx = (row * cols) + col + 2;
      const int lbl = output[idx];
      if (lbl > 0) {
        output[idx] = offset + lbl;
      }
    }
  }
}

void AddBaseOffsets(OutType &output, int cols, int num_threads, const std::vector<int> &row_starts,
                    const std::vector<int> &base) {
  tbb::parallel_for(tbb::blocked_range<int>(0, num_threads), [&](const tbb::blocked_range<int> &range) {
    for (int thread_idx = range.begin(); thread_idx < range.end(); ++thread_idx) {
      AddStripBaseOffset(thread_idx, output, cols, row_starts, base);
    }
  });
}

void MergeBoundariesUnionFind(const InType &input, const OutType &output, int rows, int cols, int num_threads,
                              const std::vector<int> &row_starts, std::vector<int> &parent) {
  for (int thread_idx = 0; thread_idx + 1 < num_threads; thread_idx++) {
    const int boundary_row = row_starts[static_cast<std::size_t>(thread_idx) + 1U];
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

std::vector<int> BuildRemapFromFirstPositions(const std::vector<int> &first_pos) {
  std::vector<int> roots;
  roots.reserve(first_pos.size());
  for (std::size_t root = 1; root < first_pos.size(); root++) {
    if (first_pos[root] != std::numeric_limits<int>::max()) {
      roots.push_back(static_cast<int>(root));
    }
  }
  std::ranges::sort(roots, [&](int a, int b) {
    return first_pos[static_cast<std::size_t>(a)] < first_pos[static_cast<std::size_t>(b)];
  });

  std::vector<int> remap(first_pos.size(), 0);
  int next_final = 1;
  for (int root : roots) {
    remap[static_cast<std::size_t>(root)] = next_final++;
  }
  return remap;
}

void CollectFirstPositionsForRows(const InType &input, const OutType &output, const std::vector<int> &parent, int cols,
                                  int row_begin, int row_end, std::vector<int> &local) {
  for (int row = row_begin; row < row_end; row++) {
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
      const int pos = (row * cols) + col;
      local[static_cast<std::size_t>(root)] = std::min(local[static_cast<std::size_t>(root)], pos);
    }
  }
}

void CollectFirstPositionsRange(const InType &input, const OutType &output, const std::vector<int> &parent, int rows,
                                int cols, int num_threads, std::vector<std::vector<int>> &thread_first, int begin_tid,
                                int end_tid) {
  for (int tid = begin_tid; tid < end_tid; ++tid) {
    auto &local = thread_first[static_cast<std::size_t>(tid)];
    const int row_begin = (tid * rows) / num_threads;
    const int row_end = ((tid + 1) * rows) / num_threads;
    CollectFirstPositionsForRows(input, output, parent, cols, row_begin, row_end, local);
  }
}

void MergeRootFirstPositions(const std::vector<std::vector<int>> &thread_first, int num_threads, int root,
                             std::vector<int> &first_pos) {
  int min_pos = std::numeric_limits<int>::max();
  for (int tid = 0; tid < num_threads; tid++) {
    min_pos = std::min(min_pos, thread_first[static_cast<std::size_t>(tid)][static_cast<std::size_t>(root)]);
  }
  first_pos[static_cast<std::size_t>(root)] = min_pos;
}

std::vector<int> ComputeFirstPositionsParallel(const InType &input, const OutType &output,
                                               const std::vector<int> &parent, int rows, int cols, int max_label,
                                               int num_threads) {
  std::vector<std::vector<int>> thread_first(
      static_cast<std::size_t>(num_threads),
      std::vector<int>(static_cast<std::size_t>(max_label) + 1U, std::numeric_limits<int>::max()));

  tbb::task_arena arena(num_threads);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<int>(0, num_threads), [&](const tbb::blocked_range<int> &range) {
      CollectFirstPositionsRange(input, output, parent, rows, cols, num_threads, thread_first, range.begin(),
                                 range.end());
    });
  });

  std::vector<int> first_pos(static_cast<std::size_t>(max_label) + 1U, std::numeric_limits<int>::max());
  tbb::parallel_for(tbb::blocked_range<int>(1, max_label + 1), [&](const tbb::blocked_range<int> &range) {
    for (int root = range.begin(); root < range.end(); ++root) {
      MergeRootFirstPositions(thread_first, num_threads, root, first_pos);
    }
  });
  return first_pos;
}

void RemapRow(const InType &input, OutType &output, const std::vector<int> &parent, const std::vector<int> &remap,
              int row, int cols) {
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
    output[idx] = remap[static_cast<std::size_t>(root)];
  }
}

void ApplyRemapInParallel(const InType &input, OutType &output, const std::vector<int> &parent,
                          const std::vector<int> &remap, int rows, int cols, int num_threads) {
  tbb::task_arena arena(num_threads);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<int>(0, rows), [&](const tbb::blocked_range<int> &range) {
      for (int row = range.begin(); row < range.end(); ++row) {
        RemapRow(input, output, parent, remap, row, cols);
      }
    });
  });
}

}  // namespace

bool GaivoronskiyMMarkingBinaryComponentsTBB::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  const int rows = input[0];
  const int cols = input[1];
  const int cells = rows * cols;

  if (cells == 0) {
    return true;
  }

  const int num_threads = NumThreadsForRows(rows);
  if (num_threads <= 1) {
    RunSequentialMarking(input, output, rows, cols);
    return true;
  }

  const std::vector<int> row_starts = MakeRowStarts(rows, num_threads);
  std::vector<int> labels_used(static_cast<std::size_t>(num_threads), 0);

  ParallelLabelStrips(input, output, cols, num_threads, row_starts, labels_used);

  const auto [base, max_global_before_merge] = BuildLabelBases(labels_used, num_threads);
  if (max_global_before_merge == 0) {
    return true;
  }

  AddBaseOffsets(output, cols, num_threads, row_starts, base);

  std::vector<int> parent(static_cast<std::size_t>(max_global_before_merge) + 1U);
  std::iota(parent.begin() + 1, parent.end(), 1);

  MergeBoundariesUnionFind(input, output, rows, cols, num_threads, row_starts, parent);
  FlattenParentForest(parent, max_global_before_merge);

  const std::vector<int> first_pos =
      ComputeFirstPositionsParallel(input, output, parent, rows, cols, max_global_before_merge, num_threads);
  const std::vector<int> remap = BuildRemapFromFirstPositions(first_pos);
  ApplyRemapInParallel(input, output, parent, remap, rows, cols, num_threads);

  return true;
}

bool GaivoronskiyMMarkingBinaryComponentsTBB::PostProcessingImpl() {
  return true;
}

}  // namespace gaivoronskiy_m_marking_binary_components
