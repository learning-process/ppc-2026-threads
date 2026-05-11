#include "zyazeva_s_matrix_mult_cannon_alg/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"

namespace zyazeva_s_matrix_mult_cannon_alg {

namespace {

using AlignedVector = std::vector<double, tbb::cache_aligned_allocator<double>>;

inline size_t BlockIndex(size_t row, size_t col, size_t grid_size) {
  return (row * grid_size) + col;
}

inline size_t BlockOffset(size_t row, size_t col, size_t grid_size, size_t block_area) {
  return BlockIndex(row, col, grid_size) * block_area;
}

size_t FindGridSize(int sz) {
  const auto max_threads = tbb::this_task_arena::max_concurrency();
  const int root = static_cast<int>(std::sqrt(max_threads));

  for (int k = root; k >= 1; --k) {
    if (sz % k == 0) {
      return static_cast<size_t>(k);
    }
  }
  return 1;
}

void CopyBlockRow(const std::vector<double> &src_mat, AlignedVector &dst_mat, size_t gi, size_t bj, size_t sz,
                  size_t bs, size_t dst_base, size_t i) {
  const size_t src = (gi * sz) + (bj * bs);
  const size_t dst = dst_base + (i * bs);
  std::copy_n(src_mat.data() + src, bs, dst_mat.data() + dst);
}

void InitializeBlock(size_t id, const std::vector<double> &m1, const std::vector<double> &m2, AlignedVector &a,
                     AlignedVector &b, size_t gs, size_t bs, size_t sz, size_t block_area) {
  const size_t bi = id / gs;
  const size_t bj = id % gs;
  const size_t base = BlockOffset(bi, bj, gs, block_area);

  for (size_t i = 0; i < bs; ++i) {
    const size_t gi = (bi * bs) + i;
    CopyBlockRow(m1, a, gi, bj, sz, bs, base, i);
    CopyBlockRow(m2, b, gi, bj, sz, bs, base, i);
  }
}

void InitializeBlocks(const std::vector<double> &m1, const std::vector<double> &m2, AlignedVector &a, AlignedVector &b,
                      size_t gs, size_t bs, size_t sz, size_t block_area, size_t total_blocks) {
  tbb::parallel_for(static_cast<size_t>(0), total_blocks,
                    [&](size_t id) { InitializeBlock(id, m1, m2, a, b, gs, bs, sz, block_area); });
}

void SetMapEntry(std::vector<size_t> &map_a, std::vector<size_t> &map_b, size_t i, size_t j, size_t gs) {
  const size_t idx = BlockIndex(i, j, gs);
  map_a[idx] = BlockIndex(i, (j + i) % gs, gs);
  map_b[idx] = BlockIndex((i + j) % gs, j, gs);
}

void InitializeMaps(std::vector<size_t> &map_a, std::vector<size_t> &map_b, size_t gs) {
  for (size_t i = 0; i < gs; ++i) {
    for (size_t j = 0; j < gs; ++j) {
      SetMapEntry(map_a, map_b, i, j, gs);
    }
  }
}

void UpdateMapEntry(std::vector<size_t> &next_a, std::vector<size_t> &next_b, const std::vector<size_t> &map_a,
                    const std::vector<size_t> &map_b, size_t i, size_t j, size_t gs) {
  const size_t idx = BlockIndex(i, j, gs);
  next_a[idx] = map_a[BlockIndex(i, (j + 1) % gs, gs)];
  next_b[idx] = map_b[BlockIndex((i + 1) % gs, j, gs)];
}

void UpdateMaps(std::vector<size_t> &map_a, std::vector<size_t> &map_b, size_t gs) {
  std::vector<size_t> next_a(map_a.size());
  std::vector<size_t> next_b(map_b.size());

  for (size_t i = 0; i < gs; ++i) {
    for (size_t j = 0; j < gs; ++j) {
      UpdateMapEntry(next_a, next_b, map_a, map_b, i, j, gs);
    }
  }

  map_a.swap(next_a);
  map_b.swap(next_b);
}

void MultiplyBlockRow(const double *a, const double *b, double *c, int bs, int i) {
  const double *a_row = a + (static_cast<size_t>(i) * bs);
  double *c_row = c + (static_cast<size_t>(i) * bs);

  for (int k = 0; k < bs; ++k) {
    const double a_val = a_row[k];
    if (a_val == 0.0) {
      continue;
    }
    const double *b_row = b + (static_cast<size_t>(k) * bs);

    for (int j = 0; j < bs; ++j) {
      c_row[j] += a_val * b_row[j];
    }
  }
}

void MultiplyBlocks(const double *a, const double *b, double *c, int block_size) {
  const int bs = block_size;
  for (int i = 0; i < bs; ++i) {
    MultiplyBlockRow(a, b, c, bs, i);
  }
}

void PerformCannonStepForBlock(size_t id, const AlignedVector &a, const AlignedVector &b, AlignedVector &c,
                               const std::vector<size_t> &map_a, const std::vector<size_t> &map_b, size_t block_area,
                               int bs) {
  const size_t a_idx = map_a[id];
  const size_t b_idx = map_b[id];

  MultiplyBlocks(a.data() + (a_idx * block_area), b.data() + (b_idx * block_area), c.data() + (id * block_area), bs);
}

void PerformCannonStep(const AlignedVector &a, const AlignedVector &b, AlignedVector &c,
                       const std::vector<size_t> &map_a, const std::vector<size_t> &map_b, size_t total_blocks,
                       size_t block_area, int bs) {
  tbb::parallel_for(static_cast<size_t>(0), total_blocks,
                    [&](size_t id) { PerformCannonStepForBlock(id, a, b, c, map_a, map_b, block_area, bs); });
}

void AssembleResultBlock(size_t id, const AlignedVector &c, std::vector<double> &result, size_t gs, size_t bs,
                         size_t sz, size_t block_area) {
  const size_t bi = id / gs;
  const size_t bj = id % gs;
  const size_t base = id * block_area;

  for (size_t i = 0; i < bs; ++i) {
    const size_t dst = ((bi * bs + i) * sz) + (bj * bs);
    const size_t src = base + (i * bs);
    std::copy_n(c.data() + src, bs, result.data() + dst);
  }
}

void AssembleResult(const AlignedVector &c, std::vector<double> &result, size_t gs, size_t bs, size_t sz,
                    size_t block_area, size_t total_blocks) {
  tbb::parallel_for(static_cast<size_t>(0), total_blocks,
                    [&](size_t id) { AssembleResultBlock(id, c, result, gs, bs, sz, block_area); });
}

}  // namespace

ZyazevaSMatrixMultCannonAlgTBB::ZyazevaSMatrixMultCannonAlgTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZyazevaSMatrixMultCannonAlgTBB::ValidationImpl() {
  const auto &input = GetInput();
  const size_t sz = std::get<0>(input);
  const auto &m1 = std::get<1>(input);
  const auto &m2 = std::get<2>(input);

  return sz > 0 && m1.size() == sz * sz && m2.size() == sz * sz;
}

bool ZyazevaSMatrixMultCannonAlgTBB::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool ZyazevaSMatrixMultCannonAlgTBB::PostProcessingImpl() {
  return true;
}

bool ZyazevaSMatrixMultCannonAlgTBB::RunImpl() {
  const int sz = static_cast<int>(std::get<0>(GetInput()));
  const auto &m1 = std::get<1>(GetInput());
  const auto &m2 = std::get<2>(GetInput());

  const size_t grid_size = FindGridSize(sz);
  const size_t bs = static_cast<size_t>(sz) / grid_size;
  const size_t gs = grid_size;
  const size_t block_area = bs * bs;
  const size_t total_blocks = gs * gs;

  AlignedVector a(total_blocks * block_area);
  AlignedVector b(total_blocks * block_area);
  AlignedVector c(total_blocks * block_area, 0.0);

  InitializeBlocks(m1, m2, a, b, gs, bs, static_cast<size_t>(sz), block_area, total_blocks);

  std::vector<size_t> map_a(total_blocks);
  std::vector<size_t> map_b(total_blocks);
  InitializeMaps(map_a, map_b, gs);

  for (size_t step = 0; step < gs; ++step) {
    PerformCannonStep(a, b, c, map_a, map_b, total_blocks, block_area, static_cast<int>(bs));
    if (step + 1 < gs) {
      UpdateMaps(map_a, map_b, gs);
    }
  }

  std::vector<double> result(static_cast<size_t>(sz) * sz);
  AssembleResult(c, result, gs, bs, static_cast<size_t>(sz), block_area, total_blocks);

  GetOutput() = std::move(result);
  return true;
}

}  // namespace zyazeva_s_matrix_mult_cannon_alg
