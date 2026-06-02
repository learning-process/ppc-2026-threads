#include "kulik_a_mat_mul_double_ccs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kulik_a_mat_mul_double_ccs {

namespace {

constexpr int kTagMatrix = 1001;

inline int ToIntCount(size_t value) {
  if (value > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("MPI count overflow");
  }
  return static_cast<int>(value);
}

inline size_t EstimateColumnCost(const CCS &a, const CCS &b, size_t j) {
  size_t cost = 0;
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    const size_t a_col = b.row[k];
    cost += a.col_ind[a_col + 1] - a.col_ind[a_col];
  }
  return cost;
}

inline std::vector<int> BuildBalancedStarts(const std::vector<size_t> &weights, int parts) {
  if (parts <= 0) {
    parts = 1;
  }
  const int n = static_cast<int>(weights.size());
  std::vector<int> starts(static_cast<size_t>(parts) + 1, 0);
  if (n == 0) {
    return starts;
  }
  std::vector<size_t> prefix(weights.size() + 1, 0);
  for (size_t i = 0; i < weights.size(); ++i) {
    prefix[i + 1] = prefix[i] + weights[i];
  }
  const size_t total = prefix.back();
  int cur = 0;
  for (int part = 1; part < parts; ++part) {
    const size_t target = total * static_cast<size_t>(part) / static_cast<size_t>(parts);
    while (cur < n && prefix[static_cast<size_t>(cur)] < target) {
      ++cur;
    }
    starts[static_cast<size_t>(part)] = cur;
  }
  starts[static_cast<size_t>(parts)] = n;
  for (int part = 1; part <= parts; ++part) {
    starts[part] = std::max(starts[part], starts[part - 1]);
  }
  return starts;
}

inline void MatMultPhase1(size_t j, const CCS &a, const CCS &b, std::vector<size_t> &was,
                          std::vector<size_t> &col_nnz) {
  size_t count = 0;
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    const size_t b_row = b.row[k];
    for (size_t zc = a.col_ind[b_row]; zc < a.col_ind[b_row + 1]; ++zc) {
      const size_t a_row = a.row[zc];
      if (was[a_row] != j) {
        was[a_row] = j;
        ++count;
      }
    }
  }
  col_nnz[j] = count;
}

inline void MatMultPhase2(size_t j, const CCS &a, const CCS &b, CCS &c, size_t stamp, std::vector<size_t> &was,
                          std::vector<double> &accum, std::vector<size_t> &rows) {
  rows.clear();
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    const double b_val = b.value[k];
    const size_t b_row = b.row[k];
    for (size_t zc = a.col_ind[b_row]; zc < a.col_ind[b_row + 1]; ++zc) {
      const size_t a_row = a.row[zc];
      accum[a_row] += a.value[zc] * b_val;
      if (was[a_row] != stamp) {
        was[a_row] = stamp;
        rows.push_back(a_row);
      }
    }
  }
  std::ranges::sort(rows);
  size_t write = c.col_ind[j];
  for (const size_t i : rows) {
    c.row[write] = i;
    c.value[write] = accum[i];
    accum[i] = 0.0;
    ++write;
  }
}

void BcastCCS(CCS &m, int root_rank = 0) {
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  MPI_Bcast(&m.m, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
  MPI_Bcast(&m.n, 1, MPI_INT, root_rank, MPI_COMM_WORLD);

  int nz = 0;
  if (world_rank == root_rank) {
    nz = static_cast<int>(m.value.size());
  }
  MPI_Bcast(&nz, 1, MPI_INT, root_rank, MPI_COMM_WORLD);

  if (world_rank != root_rank) {
    m.col_ind.resize(static_cast<size_t>(m.m) + 1);
    m.row.resize(static_cast<size_t>(nz));
    m.value.resize(static_cast<size_t>(nz));
  }

  MPI_Bcast(m.col_ind.data(), ToIntCount((static_cast<size_t>(m.m) + 1) * sizeof(size_t)), MPI_BYTE, root_rank,
            MPI_COMM_WORLD);
  MPI_Bcast(m.row.data(), ToIntCount(static_cast<size_t>(nz) * sizeof(size_t)), MPI_BYTE, root_rank, MPI_COMM_WORLD);
  MPI_Bcast(m.value.data(), nz, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
}

void SendCCS(const CCS &m, int dest) {
  MPI_Send(&m.m, 1, MPI_INT, dest, kTagMatrix, MPI_COMM_WORLD);
  MPI_Send(&m.n, 1, MPI_INT, dest, kTagMatrix, MPI_COMM_WORLD);
  const int nz = static_cast<int>(m.value.size());
  MPI_Send(&nz, 1, MPI_INT, dest, kTagMatrix, MPI_COMM_WORLD);
  MPI_Send(m.col_ind.data(), ToIntCount(m.col_ind.size() * sizeof(size_t)), MPI_BYTE, dest, kTagMatrix, MPI_COMM_WORLD);
  MPI_Send(m.row.data(), ToIntCount(m.row.size() * sizeof(size_t)), MPI_BYTE, dest, kTagMatrix, MPI_COMM_WORLD);
  MPI_Send(m.value.data(), nz, MPI_DOUBLE, dest, kTagMatrix, MPI_COMM_WORLD);
}

void RecvCCS(CCS &m, int src) {
  MPI_Status st;
  MPI_Recv(&m.m, 1, MPI_INT, src, kTagMatrix, MPI_COMM_WORLD, &st);
  MPI_Recv(&m.n, 1, MPI_INT, src, kTagMatrix, MPI_COMM_WORLD, &st);
  int nz = 0;
  MPI_Recv(&nz, 1, MPI_INT, src, kTagMatrix, MPI_COMM_WORLD, &st);
  m.col_ind.resize(static_cast<size_t>(m.m) + 1);
  m.row.resize(static_cast<size_t>(nz));
  m.value.resize(static_cast<size_t>(nz));
  MPI_Recv(m.col_ind.data(), ToIntCount(m.col_ind.size() * sizeof(size_t)), MPI_BYTE, src, kTagMatrix, MPI_COMM_WORLD,
           &st);
  MPI_Recv(m.row.data(), ToIntCount(m.row.size() * sizeof(size_t)), MPI_BYTE, src, kTagMatrix, MPI_COMM_WORLD, &st);
  MPI_Recv(m.value.data(), nz, MPI_DOUBLE, src, kTagMatrix, MPI_COMM_WORLD, &st);
}

void ScatterB(const CCS &b, CCS &b_local, const std::vector<int> &col_starts, int rank, int size) {
  if (rank == 0) {
    for (int proc = 0; proc < size; ++proc) {
      const int jstart = col_starts[static_cast<size_t>(proc)];
      const int jend = col_starts[static_cast<size_t>(proc) + 1];
      const int local_cols = jend - jstart;

      CCS tmp;
      tmp.m = local_cols;
      tmp.n = b.n;
      tmp.col_ind.resize(static_cast<size_t>(local_cols) + 1);

      const size_t nnz_start = b.col_ind[static_cast<size_t>(jstart)];
      const size_t nnz_end = b.col_ind[static_cast<size_t>(jend)];

      tmp.row.assign(b.row.begin() + static_cast<std::ptrdiff_t>(nnz_start),
                     b.row.begin() + static_cast<std::ptrdiff_t>(nnz_end));
      tmp.value.assign(b.value.begin() + static_cast<std::ptrdiff_t>(nnz_start),
                       b.value.begin() + static_cast<std::ptrdiff_t>(nnz_end));

      for (int j = 0; j <= local_cols; ++j) {
        tmp.col_ind[static_cast<size_t>(j)] =
            b.col_ind[static_cast<size_t>(jstart) + static_cast<size_t>(j)] - nnz_start;
      }

      if (proc == 0) {
        b_local = tmp;
      } else {
        SendCCS(tmp, proc);
      }
    }
  } else {
    RecvCCS(b_local, 0);
  }
}

void GatherC(CCS &c, const CCS &c_local, const std::vector<int> &col_starts, int rank, int size) {
  if (rank != 0) {
    SendCCS(c_local, 0);
    return;
  }

  std::vector<CCS> parts(static_cast<size_t>(size));
  parts[0] = c_local;

  size_t total_nnz = c_local.value.size();
  for (int proc = 1; proc < size; ++proc) {
    RecvCCS(parts[static_cast<size_t>(proc)], proc);
    total_nnz += parts[static_cast<size_t>(proc)].value.size();
  }

  c.m = col_starts.back();
  c.n = c_local.n;
  c.col_ind.assign(static_cast<size_t>(c.m) + 1, 0);
  c.row.resize(total_nnz);
  c.value.resize(total_nnz);

  size_t nnz_offset = 0;
  for (int proc = 0; proc < size; ++proc) {
    const CCS &src = parts[static_cast<size_t>(proc)];
    const int start = col_starts[static_cast<size_t>(proc)];
    const size_t cols = src.m;

    for (size_t j = 0; j <= cols; ++j) {
      c.col_ind[static_cast<size_t>(start) + j] = nnz_offset + src.col_ind[j];
    }

    std::ranges::copy(src.row, c.row.begin() + static_cast<std::ptrdiff_t>(nnz_offset));
    std::ranges::copy(src.value, c.value.begin() + static_cast<std::ptrdiff_t>(nnz_offset));
    nnz_offset += src.value.size();
  }

  c.col_ind[static_cast<size_t>(c.m)] = total_nnz;
}

}  // namespace

KulikAMatMulDoubleCcsALL::KulikAMatMulDoubleCcsALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    GetInput() = in;
    GetOutput() = CCS();
  } else {
    GetInput() = std::make_tuple(CCS(), CCS());
    GetOutput() = CCS();
  }
}

bool KulikAMatMulDoubleCcsALL::ValidationImpl() {
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    const auto &a = std::get<0>(GetInput());
    const auto &b = std::get<1>(GetInput());
    return (a.m == b.n);
  }

  return true;
}

bool KulikAMatMulDoubleCcsALL::PreProcessingImpl() {
  return true;
}

bool KulikAMatMulDoubleCcsALL::RunImpl() {
  int world_rank = 0;
  int world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  auto &a = std::get<0>(GetInput());
  auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();

  std::vector<int> col_starts;
  if (world_rank == 0) {
    std::vector<size_t> weights(static_cast<size_t>(b.m), 0);
    for (size_t j = 0; j < static_cast<size_t>(b.m); ++j) {
      weights[j] = EstimateColumnCost(a, b, j);
    }
    col_starts = BuildBalancedStarts(weights, world_size);
  } else {
    col_starts.resize(static_cast<size_t>(world_size) + 1, 0);
  }
  MPI_Bcast(col_starts.data(), world_size + 1, MPI_INT, 0, MPI_COMM_WORLD);

  BcastCCS(a);

  CCS local_b;
  ScatterB(b, local_b, col_starts, world_rank, world_size);

  const auto local_cols = static_cast<size_t>(local_b.m);

  CCS local_c;
  local_c.m = local_b.m;
  local_c.n = a.n;
  local_c.col_ind.assign(local_cols + 1, 0);

  std::vector<size_t> col_nnz(local_cols, 0);

#pragma omp parallel num_threads(std::max(1, ppc::util::GetNumThreads())) default(none) \
    shared(a, local_b, col_nnz, local_cols)
  {
    std::vector<size_t> was(static_cast<size_t>(a.n), std::numeric_limits<size_t>::max());
#pragma omp for schedule(static)
    for (size_t j = 0; j < local_cols; ++j) {
      MatMultPhase1(j, a, local_b, was, col_nnz);
    }
  }

  size_t total_nz = 0;
  for (size_t j = 0; j < local_cols; ++j) {
    local_c.col_ind[j] = total_nz;
    total_nz += col_nnz[j];
  }
  local_c.col_ind[local_cols] = total_nz;
  local_c.nz = total_nz;
  local_c.value.resize(total_nz);
  local_c.row.resize(total_nz);

#pragma omp parallel num_threads(std::max(1, ppc::util::GetNumThreads())) default(none) \
    shared(a, local_b, local_c, local_cols)
  {
    std::vector<size_t> was(static_cast<size_t>(a.n), std::numeric_limits<size_t>::max());
    std::vector<double> accum(static_cast<size_t>(a.n), 0.0);
    std::vector<size_t> rows;
#pragma omp for schedule(static)
    for (size_t j = 0; j < local_cols; ++j) {
      MatMultPhase2(j, a, local_b, local_c, local_cols + j, was, accum, rows);
    }
  }

  GatherC(c, local_c, col_starts, world_rank, world_size);

  return true;
}

bool KulikAMatMulDoubleCcsALL::PostProcessingImpl() {
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  OutType &c = GetOutput();

  int m = 0;
  int n = 0;
  int nz = 0;

  if (world_rank == 0) {
    m = static_cast<int>(c.m);
    n = static_cast<int>(c.n);
    nz = static_cast<int>(c.value.size());
  }

  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (world_rank != 0) {
    c.m = m;
    c.n = n;
    c.col_ind.resize(static_cast<size_t>(m) + 1);
    c.row.resize(static_cast<size_t>(nz));
    c.value.resize(static_cast<size_t>(nz));
  }

  MPI_Bcast(c.col_ind.data(), ToIntCount((static_cast<size_t>(m) + 1) * sizeof(size_t)), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(c.row.data(), ToIntCount(static_cast<size_t>(nz) * sizeof(size_t)), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(c.value.data(), nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return true;
}

}  // namespace kulik_a_mat_mul_double_ccs
