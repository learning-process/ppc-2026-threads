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
  starts[0] = 0;

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

inline void ProcessColumn(size_t j, const CCS &a, const CCS &b, std::vector<double> &accum,
                          std::vector<bool> &nz_elem_rows, std::vector<size_t> &nnz_rows,
                          std::vector<std::vector<double>> &local_values,
                          std::vector<std::vector<size_t>> &local_rows) {
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    const size_t ind = b.row[k];
    const double b_val = b.value[k];
    for (size_t zc = a.col_ind[ind]; zc < a.col_ind[ind + 1]; ++zc) {
      const size_t i = a.row[zc];
      const double a_val = a.value[zc];

      accum[i] += a_val * b_val;
      if (!nz_elem_rows[i]) {
        nz_elem_rows[i] = true;
        nnz_rows.push_back(i);
      }
    }
  }

  std::ranges::sort(nnz_rows);

  for (size_t i : nnz_rows) {
    if (accum[i] != 0.0) {
      local_rows[j].push_back(i);
      local_values[j].push_back(accum[i]);
    }
    accum[i] = 0.0;
    nz_elem_rows[i] = false;
  }
  nnz_rows.clear();
}

inline void CopyColumn(size_t j, CCS &c, const std::vector<std::vector<double>> &local_values,
                       const std::vector<std::vector<size_t>> &local_rows) {
  const size_t offset = c.col_ind[j];
  const size_t col_nz = local_values[j].size();
  for (size_t k = 0; k < col_nz; ++k) {
    c.value[offset + k] = local_values[j][k];
    c.row[offset + k] = local_rows[j][k];
  }
}

void ProcessColumnsRange(size_t jstart, size_t jend, const CCS &a, const CCS &b,
                         std::vector<std::vector<double>> &local_values, std::vector<std::vector<size_t>> &local_rows) {
  std::vector<double> accum(a.n, 0.0);
  std::vector<bool> nz_elem_rows(a.n, false);
  std::vector<size_t> nnz_rows;
  nnz_rows.reserve(a.n);

  for (size_t j = jstart; j < jend; ++j) {
    ProcessColumn(j, a, b, accum, nz_elem_rows, nnz_rows, local_values, local_rows);
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

  CCS local_c;
  local_c.m = local_b.m;
  local_c.n = a.n;
  local_c.col_ind.assign(static_cast<size_t>(local_c.m) + 1, 0);

  const int omp_threads = std::max(1, ppc::util::GetNumThreads());
  const int threads_count = std::max(1, std::min(omp_threads, std::max(1, static_cast<int>(local_b.m))));

  std::vector<size_t> thread_weights(static_cast<size_t>(local_b.m), 0);
  for (size_t j = 0; j < static_cast<size_t>(local_b.m); ++j) {
    thread_weights[j] = EstimateColumnCost(a, local_b, j);
  }

  std::vector<int> thread_starts = BuildBalancedStarts(thread_weights, threads_count);

  std::vector<std::vector<double>> local_values(static_cast<size_t>(local_b.m));
  std::vector<std::vector<size_t>> local_rows(static_cast<size_t>(local_b.m));

#pragma omp parallel num_threads(threads_count) default(none) \
    shared(a, local_b, local_values, local_rows, thread_starts)
  {
    const int tid = omp_get_thread_num();
    const auto jstart = static_cast<size_t>(thread_starts[static_cast<size_t>(tid)]);
    const auto jend = static_cast<size_t>(thread_starts[static_cast<size_t>(tid) + 1]);

    ProcessColumnsRange(jstart, jend, a, local_b, local_values, local_rows);
  }

  size_t total_nz = 0;
  for (size_t j = 0; j < static_cast<size_t>(local_b.m); ++j) {
    local_c.col_ind[j] = total_nz;
    total_nz += local_values[j].size();
  }
  local_c.col_ind[static_cast<size_t>(local_b.m)] = total_nz;
  local_c.value.resize(total_nz);
  local_c.row.resize(total_nz);

#pragma omp parallel for num_threads(threads_count) schedule(static) default(none) \
    shared(local_b, local_c, local_values, local_rows)
  for (size_t j = 0; j < static_cast<size_t>(local_b.m); ++j) {
    CopyColumn(j, local_c, local_values, local_rows);
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
