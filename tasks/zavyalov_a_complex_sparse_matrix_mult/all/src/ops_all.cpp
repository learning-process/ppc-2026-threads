#include "zavyalov_a_complex_sparse_matrix_mult/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {
namespace {
template <typename T>
std::vector<uint64_t> ToMPI(const std::vector<T> &v) {
  std::vector<uint64_t> res(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    res[i] = static_cast<uint64_t>(v[i]);
  }
  return res;
}

template <typename T>
std::vector<T> FromMPI(const std::vector<uint64_t> &v) {
  std::vector<T> res(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    res[i] = static_cast<T>(v[i]);
  }
  return res;
}

void BroadcastMatrix(SparseMatrix &m) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<int, 3> meta{};
  if (rank == 0) {
    meta[0] = static_cast<int>(m.height);
    meta[1] = static_cast<int>(m.width);
    meta[2] = static_cast<int>(m.val.size());
  }

  MPI_Bcast(meta.data(), 3, MPI_INT, 0, MPI_COMM_WORLD);

  m.height = static_cast<size_t>(meta[0]);
  m.width = static_cast<size_t>(meta[1]);
  auto count = static_cast<size_t>(meta[2]);

  m.row_ind.resize(count);
  m.col_ind.resize(count);
  m.val.resize(count);

  std::vector<int> tmp_rows(count);
  std::vector<int> tmp_cols(count);

  if (rank == 0) {
    for (size_t i = 0; i < count; ++i) {
      tmp_rows[i] = static_cast<int>(m.row_ind[i]);
      tmp_cols[i] = static_cast<int>(m.col_ind[i]);
    }
  }

  MPI_Bcast(tmp_rows.data(), static_cast<int>(count), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(tmp_cols.data(), static_cast<int>(count), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    for (size_t i = 0; i < count; ++i) {
      m.row_ind[i] = static_cast<size_t>(tmp_rows[i]);
      m.col_ind[i] = static_cast<size_t>(tmp_cols[i]);
    }
  }

  std::vector<double> re(count);
  std::vector<double> im(count);

  if (rank == 0) {
    for (size_t i = 0; i < count; ++i) {
      re[i] = m.val[i].re;
      im[i] = m.val[i].im;
    }
  }

  MPI_Bcast(re.data(), static_cast<int>(count), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(im.data(), static_cast<int>(count), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    for (size_t i = 0; i < count; ++i) {
      m.val[i] = Complex(re[i], im[i]);
    }
  }
}

void ScatterMatrixA(int rank, int world_size, size_t total, const SparseMatrix &ma, std::vector<int> &sendcounts,
                    std::vector<int> &displs, std::vector<size_t> &local_rows, std::vector<size_t> &local_cols,
                    std::vector<double> &local_re, std::vector<double> &local_im) {
  int blocksize = static_cast<int>(total) / world_size;
  int leftover = static_cast<int>(total) % world_size;

  int offset = 0;
  for (int proc = 0; proc < world_size; ++proc) {
    sendcounts[proc] = blocksize + (proc < leftover ? 1 : 0);
    displs[proc] = offset;
    offset += sendcounts[proc];
  }

  int local_count = sendcounts[rank];

  local_rows.resize(local_count);
  local_cols.resize(local_count);
  local_re.resize(local_count);
  local_im.resize(local_count);

  if (rank == 0) {
    for (int i = 0; i < local_count; ++i) {
      local_rows[i] = ma.row_ind[i];
      local_cols[i] = ma.col_ind[i];
      local_re[i] = ma.val[i].re;
      local_im[i] = ma.val[i].im;
    }

    for (int proc = 1; proc < world_size; ++proc) {
      int cnt = sendcounts[proc];
      int dsp = displs[proc];

      std::vector<int> rows_send(cnt);
      std::vector<int> cols_send(cnt);
      std::vector<double> re_buf(cnt);
      std::vector<double> im_buf(cnt);

      for (int i = 0; i < cnt; ++i) {
        rows_send[i] = static_cast<int>(ma.row_ind[dsp + i]);
        cols_send[i] = static_cast<int>(ma.col_ind[dsp + i]);
        re_buf[i] = ma.val[dsp + i].re;
        im_buf[i] = ma.val[dsp + i].im;
      }

      MPI_Send(rows_send.data(), cnt, MPI_INT, proc, 0, MPI_COMM_WORLD);
      MPI_Send(cols_send.data(), cnt, MPI_INT, proc, 1, MPI_COMM_WORLD);
      MPI_Send(re_buf.data(), cnt, MPI_DOUBLE, proc, 2, MPI_COMM_WORLD);
      MPI_Send(im_buf.data(), cnt, MPI_DOUBLE, proc, 3, MPI_COMM_WORLD);
    }

  } else {
    std::vector<int> rows_recv(local_count);
    std::vector<int> cols_recv(local_count);

    MPI_Recv(rows_recv.data(), local_count, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(cols_recv.data(), local_count, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(local_re.data(), local_count, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(local_im.data(), local_count, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < local_count; ++i) {
      local_rows[i] = static_cast<size_t>(rows_recv[i]);
      local_cols[i] = static_cast<size_t>(cols_recv[i]);
    }
  }
}
void GatherResult(int rank, int world_size, const std::vector<size_t> &rows, const std::vector<size_t> &cols,
                  const std::vector<double> &re_vals, const std::vector<double> &im_vals, size_t a_height,
                  size_t b_width, SparseMatrix &output) {
  int local_count = static_cast<int>(rows.size());

  std::vector<int> all_counts(world_size);
  MPI_Gather(&local_count, 1, MPI_INT, all_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displs(world_size, 0);
  int total_count = 0;

  if (rank == 0) {
    for (int i = 0; i < world_size; ++i) {
      displs[i] = total_count;
      total_count += all_counts[i];
    }
  }

  std::vector<int> all_rows(total_count);
  std::vector<int> all_cols(total_count);
  std::vector<double> all_re(total_count);
  std::vector<double> all_im(total_count);

  std::vector<int> rows_int(local_count);
  std::vector<int> cols_int(local_count);

  for (int i = 0; i < local_count; ++i) {
    rows_int[i] = static_cast<int>(rows[i]);
    cols_int[i] = static_cast<int>(cols[i]);
  }

  MPI_Gatherv(rows_int.data(), local_count, MPI_INT, all_rows.data(), all_counts.data(), displs.data(), MPI_INT, 0,
              MPI_COMM_WORLD);

  MPI_Gatherv(cols_int.data(), local_count, MPI_INT, all_cols.data(), all_counts.data(), displs.data(), MPI_INT, 0,
              MPI_COMM_WORLD);

  MPI_Gatherv(re_vals.data(), local_count, MPI_DOUBLE, all_re.data(), all_counts.data(), displs.data(), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  MPI_Gatherv(im_vals.data(), local_count, MPI_DOUBLE, all_im.data(), all_counts.data(), displs.data(), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  if (rank == 0) {
    std::map<std::pair<size_t, size_t>, Complex> mp;

    for (int i = 0; i < total_count; ++i) {
      mp[{static_cast<size_t>(all_rows[i]), static_cast<size_t>(all_cols[i])}] += Complex(all_re[i], all_im[i]);
    }

    output.height = a_height;
    output.width = b_width;

    for (auto &[k, v] : mp) {
      output.row_ind.push_back(k.first);
      output.col_ind.push_back(k.second);
      output.val.push_back(v);
    }
  }
}
}  // namespace

ZavyalovAComplSparseMatrMultALL::ZavyalovAComplSparseMatrMultALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  }
}

bool ZavyalovAComplSparseMatrMultALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());
  return matr_a.width == matr_b.height;
}

bool ZavyalovAComplSparseMatrMultALL::PreProcessingImpl() {
  return true;
}

std::map<std::pair<size_t, size_t>, Complex> ZavyalovAComplSparseMatrMultALL::ComputeLocalChunk(
    const SparseMatrix &matr_a, const SparseMatrix &matr_b, size_t start, size_t end) {
  int num_threads = ppc::util::GetNumThreads();
  std::vector<std::map<std::pair<size_t, size_t>, Complex>> local_maps(num_threads);

#pragma omp parallel for num_threads(num_threads) schedule(static) default(none) \
    shared(matr_a, matr_b, local_maps, start, end)
  for (size_t i = start; i < end; ++i) {
    int tid = omp_get_thread_num();
    size_t row_a = matr_a.row_ind[i];
    size_t col_a = matr_a.col_ind[i];
    Complex val_a = matr_a.val[i];

    for (size_t j = 0; j < matr_b.Count(); ++j) {
      if (col_a == matr_b.row_ind[j]) {
        local_maps[tid][{row_a, matr_b.col_ind[j]}] += val_a * matr_b.val[j];
      }
    }
  }

  std::map<std::pair<size_t, size_t>, Complex> result;
  for (auto &lm : local_maps) {
    for (auto &[key, value] : lm) {
      result[key] += value;
    }
  }
  return result;
}
bool ZavyalovAComplSparseMatrMultALL::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  SparseMatrix local_b = (rank == 0) ? std::get<1>(GetInput()) : SparseMatrix{};
  BroadcastMatrix(local_b);

  uint64_t total_mpi = 0;
  uint64_t a_height_mpi = 0;
  uint64_t a_width_mpi = 0;
  size_t total = 0;
  size_t a_height = 0;
  size_t a_width = 0;

  if (rank == 0) {
    const auto &ma = std::get<0>(GetInput());
    total = ma.Count();
    a_height = ma.height;
    a_width = ma.width;

    total_mpi = static_cast<uint64_t>(total);
    a_height_mpi = static_cast<uint64_t>(a_height);
    a_width_mpi = static_cast<uint64_t>(a_width);
  }

  MPI_Bcast(&total_mpi, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a_height_mpi, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a_width_mpi, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    total = static_cast<size_t>(total_mpi);
    a_height = static_cast<size_t>(a_height_mpi);
    a_width = static_cast<size_t>(a_width_mpi);
  }

  std::vector<int> sendcounts(world_size);
  std::vector<int> displs(world_size);
  std::vector<size_t> local_rows;
  std::vector<size_t> local_cols;
  std::vector<double> local_re;
  std::vector<double> local_im;

  if (rank == 0) {
    ScatterMatrixA(rank, world_size, total, std::get<0>(GetInput()), sendcounts, displs, local_rows, local_cols,
                   local_re, local_im);
  } else {
    SparseMatrix dummy;
    ScatterMatrixA(rank, world_size, total, dummy, sendcounts, displs, local_rows, local_cols, local_re, local_im);
  }

  SparseMatrix local_a;
  local_a.height = a_height;
  local_a.width = a_width;
  local_a.row_ind = local_rows;
  local_a.col_ind = local_cols;
  local_a.val.resize(local_rows.size());

  for (size_t i = 0; i < local_rows.size(); ++i) {
    local_a.val[i] = Complex(local_re[i], local_im[i]);
  }

  auto local_mp = ComputeLocalChunk(local_a, local_b, 0, local_rows.size());

  std::vector<size_t> rows;
  std::vector<size_t> cols;

  std::vector<double> re_vals;
  std::vector<double> im_vals;

  for (const auto &[key, val] : local_mp) {
    rows.push_back(key.first);
    cols.push_back(key.second);
    re_vals.push_back(val.re);
    im_vals.push_back(val.im);
  }

  SparseMatrix result;
  GatherResult(rank, world_size, rows, cols, re_vals, im_vals, a_height, local_b.width, result);

  if (rank == 0) {
    GetOutput() = std::move(result);
  }

  return true;
}

bool ZavyalovAComplSparseMatrMultALL::PostProcessingImpl() {
  return true;
}

}  // namespace zavyalov_a_compl_sparse_matr_mult
