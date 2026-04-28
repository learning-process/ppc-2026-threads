#include "zagryadskov_m_complex_spmm_ccs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <complex>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "util/include/util.hpp"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

ZagryadskovMComplexSpMMCCSALL::ZagryadskovMComplexSpMMCCSALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int world_rank = 0;
  int err_code = 0;
  err_code = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (err_code != MPI_SUCCESS) {
    throw std::runtime_error("MPI_Comm_rank failed");
  }
  if (world_rank == 0) {
    GetInput() = in;
    GetOutput() = CCS();
  } else {
    GetInput() = std::make_tuple(CCS(), CCS());
    GetOutput() = CCS();
  }
}

void ZagryadskovMComplexSpMMCCSALL::SpMMSymbolic(const CCS &a, const CCS &b, std::vector<int> &col_ptr, int jstart,
                                                 int jend) {
  std::vector<int> marker(a.m, -1);

  for (int j = jstart; j < jend; ++j) {
    int count = 0;

    for (int k = b.col_ptr[j]; k < b.col_ptr[j + 1]; ++k) {
      int b_row = b.row_ind[k];
      for (int zp = a.col_ptr[b_row]; zp < a.col_ptr[b_row + 1]; ++zp) {
        int a_row = a.row_ind[zp];
        if (marker[a_row] != j) {
          marker[a_row] = j;
          ++count;
        }
      }
    }
    col_ptr[j + 1] += count;
  }
}

void ZagryadskovMComplexSpMMCCSALL::SpMMKernel(const CCS &a, const CCS &b, CCS &c, const std::complex<double> &zero,
                                               std::vector<int> &rows, std::vector<std::complex<double>> &acc,
                                               std::vector<int> &marker, int j) {
  rows.clear();
  int write_ptr = c.col_ptr[j];

  for (int k = b.col_ptr[j]; k < b.col_ptr[j + 1]; ++k) {
    std::complex<double> tmpval = b.values[k];
    int b_row = b.row_ind[k];
    for (int zp = a.col_ptr[b_row]; zp < a.col_ptr[b_row + 1]; ++zp) {
      int a_row = a.row_ind[zp];
      acc[a_row] += tmpval * a.values[zp];
      if (marker[a_row] != j) {
        marker[a_row] = j;
        rows.push_back(a_row);
      }
    }
  }

  for (int r_idx : rows) {
    c.row_ind[write_ptr] = r_idx;
    c.values[write_ptr] = acc[r_idx];
    ++write_ptr;
    acc[r_idx] = zero;
  }
}

void ZagryadskovMComplexSpMMCCSALL::SpMMNumeric(const CCS &a, const CCS &b, CCS &c, const std::complex<double> &zero,
                                                int jstart, int jend) {
  std::vector<int> marker(a.m, -1);
  std::vector<std::complex<double>> acc(a.m, zero);
  std::vector<int> rows;

  for (int j = jstart; j < jend; ++j) {
    SpMMKernel(a, b, c, zero, rows, acc, marker, j);
  }
}

void ZagryadskovMComplexSpMMCCSALL::SpMM(const CCS &a, const CCS &b, CCS &c) {
  c.m = a.m;
  c.n = b.n;
  const int num_threads = ppc::util::GetNumThreads();

  std::complex<double> zero(0.0, 0.0);
  c.col_ptr.assign(c.n + 1, 0);

#pragma omp parallel default(none) shared(num_threads, a, b, c) num_threads(ppc::util::GetNumThreads())
  {
    int tid = omp_get_thread_num();
    int jstart = (tid * b.n) / num_threads;
    int jend = ((tid + 1) * b.n) / num_threads;
    SpMMSymbolic(a, b, c.col_ptr, jstart, jend);
  }

  for (int j = 0; j < c.n; ++j) {
    c.col_ptr[j + 1] += c.col_ptr[j];
  }
  int nnz = c.col_ptr[b.n];
  c.row_ind.resize(nnz);
  c.values.resize(nnz);
#pragma omp parallel default(none) shared(num_threads, a, b, c, zero) num_threads(ppc::util::GetNumThreads())
  {
    int tid = omp_get_thread_num();
    int jstart = (tid * b.n) / num_threads;
    int jend = ((tid + 1) * b.n) / num_threads;
    SpMMNumeric(a, b, c, zero, jstart, jend);
  }
}

void ZagryadskovMComplexSpMMCCSALL::BcastCCS(CCS &a, int rank) {
  MPI_Bcast(&a.m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a.n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int nz = 0;
  if (rank == 0) {
    nz = static_cast<int>(a.values.size());
  }
  MPI_Bcast(&nz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    a.col_ptr.resize(a.n + 1);
    a.row_ind.resize(nz);
    a.values.resize(nz);
  }

  MPI_Bcast(a.col_ptr.data(), a.n + 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(a.row_ind.data(), nz, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(a.values.data(), nz, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}

void ZagryadskovMComplexSpMMCCSALL::SendCCS(const CCS &m, int dest) {
  MPI_Send(&m.m, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
  MPI_Send(&m.n, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
  int nz = static_cast<int>(m.values.size());
  MPI_Send(&nz, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

  MPI_Send(m.col_ptr.data(), m.n + 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
  MPI_Send(m.row_ind.data(), nz, MPI_INT, dest, 0, MPI_COMM_WORLD);
  MPI_Send(m.values.data(), nz, MPI_C_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
}

void ZagryadskovMComplexSpMMCCSALL::RecvCCS(CCS &m, int src) {
  MPI_Status st;
  MPI_Recv(&m.m, 1, MPI_INT, src, 0, MPI_COMM_WORLD, &st);
  MPI_Recv(&m.n, 1, MPI_INT, src, 0, MPI_COMM_WORLD, &st);
  int nz = 0;
  MPI_Recv(&nz, 1, MPI_INT, src, 0, MPI_COMM_WORLD, &st);

  m.col_ptr.resize(m.n + 1);
  m.row_ind.resize(nz);
  m.values.resize(nz);

  MPI_Recv(m.col_ptr.data(), m.n + 1, MPI_INT, src, 0, MPI_COMM_WORLD, &st);
  MPI_Recv(m.row_ind.data(), nz, MPI_INT, src, 0, MPI_COMM_WORLD, &st);
  MPI_Recv(m.values.data(), nz, MPI_C_DOUBLE_COMPLEX, src, 0, MPI_COMM_WORLD, &st);
}

void ZagryadskovMComplexSpMMCCSALL::ScatterB(const CCS &b, CCS &b_local, const std::vector<int> &col_starts, int rank,
                                             int size) {
  if (rank == 0) {
    CCS tmp;
    for (int proc = 0; proc < size; ++proc) {
      int jstart = col_starts[proc];
      int jend = col_starts[proc + 1];

      tmp.m = b.m;
      tmp.n = jend - jstart;
      tmp.row_ind.clear();
      tmp.values.clear();
      tmp.col_ptr.clear();

      int nnz_start = b.col_ptr[jstart];
      int nnz_end = b.col_ptr[jend];
      tmp.row_ind.assign(b.row_ind.begin() + nnz_start, b.row_ind.begin() + nnz_end);
      tmp.values.assign(b.values.begin() + nnz_start, b.values.begin() + nnz_end);
      tmp.col_ptr.resize(tmp.n + 1);
      for (int j = 0; j <= tmp.n; ++j) {
        tmp.col_ptr[j] = b.col_ptr[jstart + j] - nnz_start;
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

void ZagryadskovMComplexSpMMCCSALL::GatherC(CCS &c, CCS &c_local, int rank, int size) {
  MPI_Status st;
  int local_nnz = static_cast<int>(c_local.values.size());
  int total_nnz = 0;
  int local_cols = c_local.n;
  int total_cols = 0;
  std::vector<int> tmp;
  std::vector<int> recvcounts(size);
  std::vector<int> displs(size);

  MPI_Gather(&local_nnz, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    c.m = c_local.m;
    for (int i = 0; i < size; ++i) {
      displs[i] = total_nnz;
      total_nnz += recvcounts[i];
    }
    c.row_ind.resize(total_nnz);
    c.values.resize(total_nnz);
  }

  MPI_Gatherv(c_local.row_ind.data(), local_nnz, MPI_INT, c.row_ind.data(), recvcounts.data(), displs.data(), MPI_INT,
              0, MPI_COMM_WORLD);
  MPI_Gatherv(c_local.values.data(), local_nnz, MPI_C_DOUBLE_COMPLEX, c.values.data(), recvcounts.data(), displs.data(),
              MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

  MPI_Gather(&local_cols, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (int i = 0; i < size; ++i) {
      displs[i] = total_cols + 1;
      total_cols += recvcounts[i];
      recvcounts[i] += 1;
    }
    c.n = total_cols;
    c.col_ptr.resize(total_cols + 1);
  }

  if (rank != 0) {
    MPI_Send(c_local.col_ptr.data(), c_local.n + 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }

  if (rank == 0) {
    std::ranges::copy(c_local.col_ptr, c.col_ptr.begin());

    int nz_offset = c_local.col_ptr.back();
    int col_offset = c_local.n;
    for (int proc = 1; proc < size; ++proc) {
      tmp.resize(recvcounts[proc]);
      MPI_Recv(tmp.data(), recvcounts[proc], MPI_INT, proc, 0, MPI_COMM_WORLD, &st);

      for (int j = 1; j < recvcounts[proc]; ++j) {
        c.col_ptr[col_offset + j] = nz_offset + tmp[j];
      }

      nz_offset += tmp.back();
      col_offset += recvcounts[proc] - 1;
      tmp.clear();
    }
  }
}

bool ZagryadskovMComplexSpMMCCSALL::ValidationImpl() {
  bool res = false;
  int world_rank = 0;
  int err_code = 0;
  err_code = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (err_code != MPI_SUCCESS) {
    throw std::runtime_error("MPI_Comm_rank failed");
  }
  if (world_rank == 0) {
    const CCS &a = std::get<0>(GetInput());
    const CCS &b = std::get<1>(GetInput());
    res = a.n == b.m;
  } else {
    res = true;
  }
  return res;
}

bool ZagryadskovMComplexSpMMCCSALL::PreProcessingImpl() {
  return true;
}

bool ZagryadskovMComplexSpMMCCSALL::RunImpl() {
  int world_rank = 0;
  int world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  CCS &a = std::get<0>(GetInput());
  CCS &b = std::get<1>(GetInput());
  CCS &c = GetOutput();

  CCS local_b;
  CCS local_c;
  std::vector<int> col_starts;
  if (world_rank == 0) {
    col_starts.resize(world_size + 1);
    for (int proc = 0; proc <= world_size; ++proc) {
      col_starts[proc] = (proc * b.n) / world_size;
    }
  }

  BcastCCS(a, world_rank);
  ScatterB(b, local_b, col_starts, world_rank, world_size);

  ZagryadskovMComplexSpMMCCSALL::SpMM(a, local_b, local_c);

  GatherC(c, local_c, world_rank, world_size);

  return true;
}

bool ZagryadskovMComplexSpMMCCSALL::PostProcessingImpl() {
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int m = 0;
  int n = 0;
  int nz = 0;
  CCS &c = GetOutput();
  if (world_rank == 0) {
    m = c.m;
    n = c.n;
    nz = static_cast<int>(c.values.size());
  }
  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (world_rank != 0) {
    c.m = m;
    c.n = n;
    c.col_ptr.resize(n + 1);
    c.row_ind.resize(nz);
    c.values.resize(nz);
  }
  MPI_Bcast(c.col_ptr.data(), n + 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(c.row_ind.data(), nz, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(c.values.data(), nz, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

  return true;
}

}  // namespace zagryadskov_m_complex_spmm_ccs
