#include "potashnik_m_matrix_mult_complex/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <thread>
#include <utility>
#include <vector>

#include "potashnik_m_matrix_mult_complex/common/include/common.hpp"

namespace potashnik_m_matrix_mult_complex {

namespace {

using Key = std::pair<size_t, size_t>;
using LocalMap = std::map<Key, Complex>;

void BroadcastMatrix(CCSMatrix &matrix) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<int, 3> meta{};
  if (rank == 0) {
    meta[0] = static_cast<int>(matrix.height);
    meta[1] = static_cast<int>(matrix.width);
    meta[2] = static_cast<int>(matrix.val.size());
  }
  MPI_Bcast(meta.data(), 3, MPI_INT, 0, MPI_COMM_WORLD);
  matrix.height = static_cast<size_t>(meta[0]);
  matrix.width = static_cast<size_t>(meta[1]);
  auto count = static_cast<size_t>(meta[2]);

  matrix.val.resize(count);
  matrix.row_ind.resize(count);
  matrix.col_ptr.resize(count);

  std::vector<int> tmp_rows(count);
  std::vector<int> tmp_cols(count);
  std::vector<double> re(count);
  std::vector<double> im(count);

  if (rank == 0) {
    for (size_t i = 0; i < count; ++i) {
      tmp_rows[i] = static_cast<int>(matrix.row_ind[i]);
      tmp_cols[i] = static_cast<int>(matrix.col_ptr[i]);
      re[i] = matrix.val[i].real;
      im[i] = matrix.val[i].imaginary;
    }
  }

  MPI_Bcast(tmp_rows.data(), static_cast<int>(count), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(tmp_cols.data(), static_cast<int>(count), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(re.data(), static_cast<int>(count), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(im.data(), static_cast<int>(count), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    for (size_t i = 0; i < count; ++i) {
      matrix.row_ind[i] = static_cast<size_t>(tmp_rows[i]);
      matrix.col_ptr[i] = static_cast<size_t>(tmp_cols[i]);
      matrix.val[i] = Complex(re[i], im[i]);
    }
  }
}

void ScatterMatrixLeft(int rank, int world_size, size_t total, const CCSMatrix &matrix, std::vector<int> &sendcounts,
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
      local_rows[i] = matrix.row_ind[i];
      local_cols[i] = matrix.col_ptr[i];
      local_re[i] = matrix.val[i].real;
      local_im[i] = matrix.val[i].imaginary;
    }
    for (int proc = 1; proc < world_size; ++proc) {
      int cnt = sendcounts[proc];
      int dsp = displs[proc];
      std::vector<int> rows_send(cnt);
      std::vector<int> cols_send(cnt);
      std::vector<double> re_buf(cnt);
      std::vector<double> im_buf(cnt);
      for (int i = 0; i < cnt; ++i) {
        rows_send[i] = static_cast<int>(matrix.row_ind[dsp + i]);
        cols_send[i] = static_cast<int>(matrix.col_ptr[dsp + i]);
        re_buf[i] = matrix.val[dsp + i].real;
        im_buf[i] = matrix.val[dsp + i].imaginary;
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
                  const std::vector<double> &re_vals, const std::vector<double> &im_vals, size_t height_left,
                  size_t width_right, CCSMatrix &output) {
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
    std::map<Key, Complex> buffer;
    for (int i = 0; i < total_count; ++i) {
      buffer[{static_cast<size_t>(all_rows[i]), static_cast<size_t>(all_cols[i])}] += Complex(all_re[i], all_im[i]);
    }
    output.height = height_left;
    output.width = width_right;
    output.val.reserve(buffer.size());
    output.row_ind.reserve(buffer.size());
    output.col_ptr.reserve(buffer.size());
    for (const auto &[key, value] : buffer) {
      output.row_ind.push_back(key.first);
      output.col_ptr.push_back(key.second);
      output.val.push_back(value);
    }
  }
}

void ProcessChunk(size_t begin, size_t end, const CCSMatrix &matrix_right, const std::vector<Complex> &val_left,
                  const std::vector<size_t> &row_ind_left, const std::vector<size_t> &col_ptr_left,
                  LocalMap &local_buffer) {
  const auto &val_right = matrix_right.val;
  const auto &row_ind_right = matrix_right.row_ind;
  const auto &col_ptr_right = matrix_right.col_ptr;

  for (size_t i = begin; i < end; ++i) {
    size_t row_left = row_ind_left[i];
    size_t col_left = col_ptr_left[i];
    Complex left_val = val_left[i];

    for (size_t j = 0; j < matrix_right.Count(); ++j) {
      size_t row_right = row_ind_right[j];
      size_t col_right = col_ptr_right[j];
      Complex right_val = val_right[j];

      if (col_left == row_right) {
        local_buffer[{row_left, col_right}] += left_val * right_val;
      }
    }
  }
}

LocalMap ComputeLocalChunk(const CCSMatrix &matrix_left, const CCSMatrix &matrix_right) {
  const auto &val_left = matrix_left.val;
  const auto &row_ind_left = matrix_left.row_ind;
  const auto &col_ptr_left = matrix_left.col_ptr;

  size_t left_count = matrix_left.Count();
  size_t num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  std::vector<LocalMap> local_buffers(num_threads);
  std::vector<std::thread> threads(num_threads);
  size_t chunk = (left_count + num_threads - 1) / num_threads;

  for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t begin = thread_idx * chunk;
    size_t end = std::min(begin + chunk, left_count);
    threads[thread_idx] = std::thread([&, thread_idx, begin, end]() {
      ProcessChunk(begin, end, matrix_right, val_left, row_ind_left, col_ptr_left, local_buffers[thread_idx]);
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  LocalMap result;
  for (const auto &local : local_buffers) {
    for (const auto &[key, value] : local) {
      result[key] += value;
    }
  }
  return result;
}

}  // namespace

PotashnikMMatrixMultComplexALL::PotashnikMMatrixMultComplexALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  }
}

bool PotashnikMMatrixMultComplexALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }
  const auto &matrix_left = std::get<0>(GetInput());
  const auto &matrix_right = std::get<1>(GetInput());
  return matrix_left.width == matrix_right.height;
}

bool PotashnikMMatrixMultComplexALL::PreProcessingImpl() {
  return true;
}

bool PotashnikMMatrixMultComplexALL::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  CCSMatrix matrix_right = (rank == 0) ? std::get<1>(GetInput()) : CCSMatrix{};
  BroadcastMatrix(matrix_right);

  std::array<uint64_t, 3> meta{};
  if (rank == 0) {
    const auto &ml = std::get<0>(GetInput());
    meta[0] = static_cast<uint64_t>(ml.Count());
    meta[1] = static_cast<uint64_t>(ml.height);
    meta[2] = static_cast<uint64_t>(ml.width);
  }
  MPI_Bcast(meta.data(), 3, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  auto total = static_cast<size_t>(meta[0]);
  auto height_left = static_cast<size_t>(meta[1]);

  std::vector<int> sendcounts(world_size);
  std::vector<int> displs(world_size);
  std::vector<size_t> local_rows;
  std::vector<size_t> local_cols;
  std::vector<double> local_re;
  std::vector<double> local_im;

  const CCSMatrix empty{};
  const CCSMatrix &matrix_left_ref = (rank == 0) ? std::get<0>(GetInput()) : empty;
  ScatterMatrixLeft(rank, world_size, total, matrix_left_ref, sendcounts, displs, local_rows, local_cols, local_re,
                    local_im);

  CCSMatrix local_left;
  local_left.height = height_left;
  local_left.width = matrix_right.height;
  local_left.row_ind = local_rows;
  local_left.col_ptr = local_cols;
  local_left.val.resize(local_rows.size());
  for (size_t i = 0; i < local_rows.size(); ++i) {
    local_left.val[i] = Complex(local_re[i], local_im[i]);
  }

  auto local_result = ComputeLocalChunk(local_left, matrix_right);

  std::vector<size_t> res_rows;
  std::vector<size_t> res_cols;
  std::vector<double> res_re;
  std::vector<double> res_im;
  for (const auto &[key, value] : local_result) {
    res_rows.push_back(key.first);
    res_cols.push_back(key.second);
    res_re.push_back(value.real);
    res_im.push_back(value.imaginary);
  }

  CCSMatrix output;
  GatherResult(rank, world_size, res_rows, res_cols, res_re, res_im, height_left, matrix_right.width, output);

  if (rank == 0) {
    GetOutput() = std::move(output);
  }
  return true;
}

bool PotashnikMMatrixMultComplexALL::PostProcessingImpl() {
  return true;
}

}  // namespace potashnik_m_matrix_mult_complex
