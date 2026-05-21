#include "cheremkhin_a_matr_mult_cannon_alg/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "cheremkhin_a_matr_mult_cannon_alg/common/include/common.hpp"
#include "util/include/util.hpp"

namespace cheremkhin_a_matr_mult_cannon_alg {

namespace {

inline std::size_t Idx(std::size_t n, std::size_t r, std::size_t c) {
  return (r * n) + c;
}

std::size_t CeilDiv(std::size_t a, std::size_t b) {
  return (a + b - 1) / b;
}

int ChooseVirtualGridSize(int world_size) {
  if (world_size <= 1) {
    return 1;
  }

  int grid_dim = static_cast<int>(std::sqrt(static_cast<double>(world_size)));
  while ((grid_dim * grid_dim) < world_size) {
    ++grid_dim;
  }
  return grid_dim;
}

int MakeVirtualRank(int row, int col, int grid_dim) {
  return (row * grid_dim) + col;
}

int GetOwnerRank(int virtual_rank, int world_size) {
  return virtual_rank % world_size;
}

std::vector<int> GetOwnedVirtualRanks(int world_rank, int world_size, int grid_dim) {
  std::vector<int> owned_ranks;
  const int virtual_size = grid_dim * grid_dim;
  for (int virtual_rank = world_rank; virtual_rank < virtual_size; virtual_rank += world_size) {
    owned_ranks.push_back(virtual_rank);
  }
  return owned_ranks;
}

void CopyGlobalToPadded(const std::vector<double> &src, std::vector<double> &dst, std::size_t src_n,
                        std::size_t dst_n) {
  const auto src_n64 = static_cast<std::int64_t>(src_n);
#pragma omp parallel for default(none) schedule(static) shared(src, dst, src_n, dst_n, src_n64)
  for (std::int64_t i = 0; i < src_n64; ++i) {
    for (std::size_t j = 0; j < src_n; ++j) {
      dst[Idx(dst_n, static_cast<std::size_t>(i), j)] = src[Idx(src_n, static_cast<std::size_t>(i), j)];
    }
  }
}

void CopyPaddedToGlobal(const std::vector<double> &src, std::vector<double> &dst, std::size_t src_n,
                        std::size_t dst_n) {
  const auto dst_n64 = static_cast<std::int64_t>(dst_n);
#pragma omp parallel for default(none) schedule(static) shared(src, dst, src_n, dst_n, dst_n64)
  for (std::int64_t i = 0; i < dst_n64; ++i) {
    for (std::size_t j = 0; j < dst_n; ++j) {
      dst[Idx(dst_n, static_cast<std::size_t>(i), j)] = src[Idx(src_n, static_cast<std::size_t>(i), j)];
    }
  }
}

void ExtractLocalBlock(const std::vector<double> &src, std::vector<double> &block, std::size_t global_n,
                       std::size_t block_n, int block_row, int block_col) {
  const std::size_t row0 = static_cast<std::size_t>(block_row) * block_n;
  const std::size_t col0 = static_cast<std::size_t>(block_col) * block_n;
  const auto block_n64 = static_cast<std::int64_t>(block_n);
#pragma omp parallel for default(none) schedule(static) shared(src, block, global_n, block_n, row0, col0, block_n64)
  for (std::int64_t i = 0; i < block_n64; ++i) {
    const std::size_t src_row = (row0 + static_cast<std::size_t>(i)) * global_n;
    const std::size_t dst_row = static_cast<std::size_t>(i) * block_n;
    for (std::size_t j = 0; j < block_n; ++j) {
      block[dst_row + j] = src[src_row + col0 + j];
    }
  }
}

void InsertLocalBlock(const std::vector<double> &block, std::vector<double> &dst, std::size_t global_n,
                      std::size_t block_n, int block_row, int block_col) {
  const std::size_t row0 = static_cast<std::size_t>(block_row) * block_n;
  const std::size_t col0 = static_cast<std::size_t>(block_col) * block_n;
  const auto block_n64 = static_cast<std::int64_t>(block_n);
#pragma omp parallel for default(none) schedule(static) shared(block, dst, global_n, block_n, row0, col0, block_n64)
  for (std::int64_t i = 0; i < block_n64; ++i) {
    const std::size_t src_row = static_cast<std::size_t>(i) * block_n;
    const std::size_t dst_row = (row0 + static_cast<std::size_t>(i)) * global_n;
    for (std::size_t j = 0; j < block_n; ++j) {
      dst[dst_row + col0 + j] = block[src_row + j];
    }
  }
}

void MulAddLocal(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                 std::size_t block_n) {
  const auto block_n64 = static_cast<std::int64_t>(block_n);

#pragma omp parallel for default(none) schedule(static) shared(a, b, c, block_n, block_n64)
  for (std::int64_t ii = 0; ii < block_n64; ++ii) {
    const auto row = static_cast<std::size_t>(ii);
    const std::size_t a_row = row * block_n;
    const std::size_t c_row = row * block_n;
    double *c_block = c.data() + c_row;
    for (std::size_t kk = 0; kk < block_n; ++kk) {
      const double aik = a[a_row + kk];
      const double *b_block = b.data() + (kk * block_n);
      for (std::int64_t jj = 0; jj < block_n64; ++jj) {
        c_block[jj] += aik * b_block[jj];
      }
    }
  }
}

struct LocalCell {
  int virtual_rank = 0;
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;
};

int GetRow(int virtual_rank, int grid_dim) {
  return virtual_rank / grid_dim;
}

int GetCol(int virtual_rank, int grid_dim) {
  return virtual_rank % grid_dim;
}

struct ShiftTargets {
  int source_rank = 0;
  int dest_rank = 0;
  int source_owner = 0;
  int dest_owner = 0;
};

ShiftTargets ComputeShiftTargets(int virtual_rank, const std::vector<int> &owner_by_rank, int grid_dim,
                                 bool horizontal_shift) {
  const int row = GetRow(virtual_rank, grid_dim);
  const int col = GetCol(virtual_rank, grid_dim);

  ShiftTargets targets;
  targets.source_rank = horizontal_shift ? MakeVirtualRank(row, (col + 1) % grid_dim, grid_dim)
                                         : MakeVirtualRank((row + 1) % grid_dim, col, grid_dim);
  targets.dest_rank = horizontal_shift ? MakeVirtualRank(row, (col + grid_dim - 1) % grid_dim, grid_dim)
                                       : MakeVirtualRank((row + grid_dim - 1) % grid_dim, col, grid_dim);
  targets.source_owner = owner_by_rank[static_cast<std::size_t>(targets.source_rank)];
  targets.dest_owner = owner_by_rank[static_cast<std::size_t>(targets.dest_rank)];
  return targets;
}

void ExchangePhase(const std::vector<std::vector<double>> &current_buffers,
                   std::vector<std::vector<double>> &next_buffers, const std::vector<int> &virtual_ranks,
                   const std::vector<int> &owner_by_rank, const std::vector<int> &local_index_by_rank, int grid_dim,
                   int world_rank, int tag_base, bool horizontal_shift) {
  std::size_t recv_count = 0;
  std::size_t send_count = 0;
  for (int virtual_rank : virtual_ranks) {
    const auto targets = ComputeShiftTargets(virtual_rank, owner_by_rank, grid_dim, horizontal_shift);
    recv_count += (targets.source_owner != world_rank) ? 1U : 0U;
    send_count += (targets.dest_owner != world_rank) ? 1U : 0U;
  }

  std::vector<MPI_Request> recv_requests(recv_count, MPI_REQUEST_NULL);
  std::vector<MPI_Request> send_requests(send_count, MPI_REQUEST_NULL);
  std::size_t recv_idx = 0;
  std::size_t send_idx = 0;

  for (std::size_t idx = 0; idx < virtual_ranks.size(); ++idx) {
    const int virtual_rank = virtual_ranks[idx];
    const auto targets = ComputeShiftTargets(virtual_rank, owner_by_rank, grid_dim, horizontal_shift);

    if (targets.source_owner == world_rank) {
      const int local_source_idx = local_index_by_rank[static_cast<std::size_t>(targets.source_rank)];
      next_buffers[idx] = current_buffers[static_cast<std::size_t>(local_source_idx)];
    } else {
      MPI_Irecv(next_buffers[idx].data(), static_cast<int>(next_buffers[idx].size()), MPI_DOUBLE, targets.source_owner,
                tag_base + virtual_rank, MPI_COMM_WORLD, &recv_requests[recv_idx]);
      ++recv_idx;
    }

    if (targets.dest_owner != world_rank) {
      MPI_Isend(current_buffers[idx].data(), static_cast<int>(current_buffers[idx].size()), MPI_DOUBLE,
                targets.dest_owner, tag_base + targets.dest_rank, MPI_COMM_WORLD, &send_requests[send_idx]);
      ++send_idx;
    }
  }

  if (!recv_requests.empty()) {
    MPI_Waitall(static_cast<int>(recv_requests.size()), recv_requests.data(), MPI_STATUSES_IGNORE);
  }
  if (!send_requests.empty()) {
    MPI_Waitall(static_cast<int>(send_requests.size()), send_requests.data(), MPI_STATUSES_IGNORE);
  }
}

void DistributeInitiallyAlignedBlocks(const std::vector<double> &a_global, const std::vector<double> &b_global,
                                      std::vector<LocalCell> &local_cells, const std::vector<int> &local_index_by_rank,
                                      std::size_t global_n, std::size_t block_n, int grid_dim, int world_rank,
                                      int world_size) {
  constexpr int kTagA = 1000;
  constexpr int kTagB = 2000;

  if (world_rank == 0) {
    for (int row = 0; row < grid_dim; ++row) {
      for (int col = 0; col < grid_dim; ++col) {
        const int virtual_rank = MakeVirtualRank(row, col, grid_dim);
        const int owner_rank = GetOwnerRank(virtual_rank, world_size);
        const int a_col = (row + col) % grid_dim;
        const int b_row = (row + col) % grid_dim;

        std::vector<double> a_block(block_n * block_n, 0.0);
        std::vector<double> b_block(block_n * block_n, 0.0);
        ExtractLocalBlock(a_global, a_block, global_n, block_n, row, a_col);
        ExtractLocalBlock(b_global, b_block, global_n, block_n, b_row, col);

        if (owner_rank == 0) {
          const int local_idx = local_index_by_rank[static_cast<std::size_t>(virtual_rank)];
          local_cells[static_cast<std::size_t>(local_idx)].a = std::move(a_block);
          local_cells[static_cast<std::size_t>(local_idx)].b = std::move(b_block);
        } else {
          MPI_Send(a_block.data(), static_cast<int>(a_block.size()), MPI_DOUBLE, owner_rank, kTagA + virtual_rank,
                   MPI_COMM_WORLD);
          MPI_Send(b_block.data(), static_cast<int>(b_block.size()), MPI_DOUBLE, owner_rank, kTagB + virtual_rank,
                   MPI_COMM_WORLD);
        }
      }
    }
  } else {
    for (auto &cell : local_cells) {
      MPI_Recv(cell.a.data(), static_cast<int>(cell.a.size()), MPI_DOUBLE, 0, kTagA + cell.virtual_rank, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(cell.b.data(), static_cast<int>(cell.b.size()), MPI_DOUBLE, 0, kTagB + cell.virtual_rank, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }
}

void ShiftBlocksCannon(std::vector<LocalCell> &local_cells, const std::vector<int> &owner_by_rank, std::size_t block_n,
                       int grid_dim, int world_rank) {
  constexpr int kShiftATagBase = 3000;
  constexpr int kShiftBTagBase = 5000;

  std::vector<int> virtual_ranks(local_cells.size(), 0);
  std::vector<int> local_index_by_rank(static_cast<std::size_t>(grid_dim * grid_dim), -1);
  std::vector<std::vector<double>> current_a(local_cells.size());
  std::vector<std::vector<double>> current_b(local_cells.size());
  std::vector<std::vector<double>> next_a(local_cells.size(), std::vector<double>(block_n * block_n, 0.0));
  std::vector<std::vector<double>> next_b(local_cells.size(), std::vector<double>(block_n * block_n, 0.0));

  for (std::size_t idx = 0; idx < local_cells.size(); ++idx) {
    virtual_ranks[idx] = local_cells[idx].virtual_rank;
    local_index_by_rank[static_cast<std::size_t>(local_cells[idx].virtual_rank)] = static_cast<int>(idx);
    current_a[idx] = local_cells[idx].a;
    current_b[idx] = local_cells[idx].b;
  }

  ExchangePhase(current_a, next_a, virtual_ranks, owner_by_rank, local_index_by_rank, grid_dim, world_rank,
                kShiftATagBase, true);
  ExchangePhase(current_b, next_b, virtual_ranks, owner_by_rank, local_index_by_rank, grid_dim, world_rank,
                kShiftBTagBase, false);

  for (std::size_t idx = 0; idx < local_cells.size(); ++idx) {
    local_cells[idx].a = std::move(next_a[idx]);
    local_cells[idx].b = std::move(next_b[idx]);
  }
}

void GatherResultBlocks(const std::vector<LocalCell> &local_cells, std::vector<double> &global_matrix,
                        const std::vector<int> &local_index_by_rank, std::size_t global_n, std::size_t block_n,
                        int grid_dim, int world_rank, int world_size) {
  constexpr int kTagC = 7000;

  if (world_rank == 0) {
    for (int virtual_rank = 0; virtual_rank < grid_dim * grid_dim; ++virtual_rank) {
      const int row = virtual_rank / grid_dim;
      const int col = virtual_rank % grid_dim;
      const int owner_rank = GetOwnerRank(virtual_rank, world_size);

      if (owner_rank == 0) {
        const int local_idx = local_index_by_rank[static_cast<std::size_t>(virtual_rank)];
        InsertLocalBlock(local_cells[static_cast<std::size_t>(local_idx)].c, global_matrix, global_n, block_n, row,
                         col);
      } else {
        std::vector<double> block(block_n * block_n, 0.0);
        MPI_Recv(block.data(), static_cast<int>(block.size()), MPI_DOUBLE, owner_rank, kTagC + virtual_rank,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        InsertLocalBlock(block, global_matrix, global_n, block_n, row, col);
      }
    }
  } else {
    for (const auto &cell : local_cells) {
      MPI_Send(cell.c.data(), static_cast<int>(cell.c.size()), MPI_DOUBLE, 0, kTagC + cell.virtual_rank,
               MPI_COMM_WORLD);
    }
  }
}

}  // namespace

CheremkhinAMatrMultCannonAlgALL::CheremkhinAMatrMultCannonAlgALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool CheremkhinAMatrMultCannonAlgALL::ValidationImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a = std::get<1>(GetInput());
  const auto &b = std::get<2>(GetInput());
  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool CheremkhinAMatrMultCannonAlgALL::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool CheremkhinAMatrMultCannonAlgALL::RunImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a_in = std::get<1>(GetInput());
  const auto &b_in = std::get<2>(GetInput());
  const int requested_threads = ppc::util::GetNumThreads();
  int world_rank = 0;
  int world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  omp_set_num_threads(requested_threads);

  const int q = ChooseVirtualGridSize(world_size);
  const int virtual_size = q * q;
  const std::size_t block_n = CeilDiv(n, static_cast<std::size_t>(q));
  const std::size_t padded_n = block_n * static_cast<std::size_t>(q);

  std::vector<int> owner_by_rank(static_cast<std::size_t>(virtual_size), 0);
  for (int virtual_rank = 0; virtual_rank < virtual_size; ++virtual_rank) {
    owner_by_rank[static_cast<std::size_t>(virtual_rank)] = GetOwnerRank(virtual_rank, world_size);
  }

  const std::vector<int> owned_virtual_ranks = GetOwnedVirtualRanks(world_rank, world_size, q);
  std::vector<int> local_index_by_rank(static_cast<std::size_t>(virtual_size), -1);
  std::vector<LocalCell> local_cells;
  local_cells.reserve(owned_virtual_ranks.size());
  for (std::size_t idx = 0; idx < owned_virtual_ranks.size(); ++idx) {
    const int virtual_rank = owned_virtual_ranks[idx];
    local_index_by_rank[static_cast<std::size_t>(virtual_rank)] = static_cast<int>(idx);
    LocalCell cell;
    cell.virtual_rank = virtual_rank;
    cell.a.assign(block_n * block_n, 0.0);
    cell.b.assign(block_n * block_n, 0.0);
    cell.c.assign(block_n * block_n, 0.0);
    local_cells.push_back(std::move(cell));
  }

  std::vector<double> a_padded;
  std::vector<double> b_padded;
  if (world_rank == 0) {
    a_padded.assign(padded_n * padded_n, 0.0);
    b_padded.assign(padded_n * padded_n, 0.0);
    CopyGlobalToPadded(a_in, a_padded, n, padded_n);
    CopyGlobalToPadded(b_in, b_padded, n, padded_n);
  }

  DistributeInitiallyAlignedBlocks(a_padded, b_padded, local_cells, local_index_by_rank, padded_n, block_n, q,
                                   world_rank, world_size);

  for (int step = 0; step < q; ++step) {
    for (auto &cell : local_cells) {
      MulAddLocal(cell.a, cell.b, cell.c, block_n);
    }
    if (step + 1 < q) {
      ShiftBlocksCannon(local_cells, owner_by_rank, block_n, q, world_rank);
    }
  }

  std::vector<double> c_padded(padded_n * padded_n, 0.0);
  GatherResultBlocks(local_cells, c_padded, local_index_by_rank, padded_n, block_n, q, world_rank, world_size);
  MPI_Bcast(c_padded.data(), static_cast<int>(c_padded.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::vector<double> out(n * n, 0.0);
  CopyPaddedToGlobal(c_padded, out, padded_n, n);

  GetOutput() = std::move(out);
  return true;
}

bool CheremkhinAMatrMultCannonAlgALL::PostProcessingImpl() {
  return true;
}

}  // namespace  cheremkhin_a_matr_mult_cannon_alg
