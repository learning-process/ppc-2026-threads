#include "tabalaev_a_matrix_mul_strassen/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stack>
#include <utility>
#include <vector>

#include "tabalaev_a_matrix_mul_strassen/common/include/common.hpp"

namespace tabalaev_a_matrix_mul_strassen {

static constexpr size_t kBaseCaseSize = 128;

TabalaevAMatrixMulStrassenALL::TabalaevAMatrixMulStrassenALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool TabalaevAMatrixMulStrassenALL::ValidationImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int hw_concurrency = omp_get_num_procs();
  int threads_per_process = std::max(1, hw_concurrency / size);
  omp_set_num_threads(threads_per_process);

  int is_valid = 0;
  if (rank == 0) {
    const auto &in = GetInput();
    bool valid = in.a_rows > 0 && in.a_cols_b_rows > 0 && in.b_cols > 0 &&
                 in.a.size() == static_cast<size_t>(in.a_rows * in.a_cols_b_rows) &&
                 in.b.size() == static_cast<size_t>(in.a_cols_b_rows * in.b_cols);
    is_valid = valid ? 1 : 0;
  }
  MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return is_valid == 1;
}

bool TabalaevAMatrixMulStrassenALL::PreProcessingImpl() {
  GetOutput() = {};

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    const auto &in = GetInput();

    a_rows_ = in.a_rows;
    a_cols_b_rows_ = in.a_cols_b_rows;
    b_cols_ = in.b_cols;

    size_t max_dim = std::max({a_rows_, a_cols_b_rows_, b_cols_});
    padded_n_ = 1;
    while (padded_n_ < max_dim) {
      padded_n_ *= 2;
    }

    padded_a_.assign(padded_n_ * padded_n_, 0.0);
    padded_b_.assign(padded_n_ * padded_n_, 0.0);

    auto &padded_a = padded_a_;
    auto &padded_b = padded_b_;
    size_t a_rows = a_rows_;
    size_t a_cols_b_rows = a_cols_b_rows_;
    size_t b_cols = b_cols_;
    size_t padded_n = padded_n_;

#pragma omp parallel default(none) shared(in, padded_a, padded_b, a_rows, a_cols_b_rows, b_cols, padded_n)
    {
#pragma omp for nowait
      for (size_t i = 0; i < a_rows; ++i) {
        for (size_t j = 0; j < a_cols_b_rows; ++j) {
          padded_a[(i * padded_n) + j] = in.a[(i * a_cols_b_rows) + j];
        }
      }
#pragma omp for
      for (size_t i = 0; i < a_cols_b_rows; ++i) {
        for (size_t j = 0; j < b_cols; ++j) {
          padded_b[(i * padded_n) + j] = in.b[(i * b_cols) + j];
        }
      }
    }
  }

  return true;
}

bool TabalaevAMatrixMulStrassenALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int procs = omp_get_num_procs();
  int threads_per_process = std::max(1, procs / size);
  omp_set_num_threads(threads_per_process);

  std::array<uint64_t, 3> dims = {0, 0, 0};
  if (rank == 0) {
    dims[0] = static_cast<uint64_t>(a_rows_);
    dims[1] = static_cast<uint64_t>(b_cols_);
    dims[2] = static_cast<uint64_t>(padded_n_);
  }
  MPI_Bcast(dims.data(), 3, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    a_rows_ = static_cast<size_t>(dims[0]);
    b_cols_ = static_cast<size_t>(dims[1]);
    padded_n_ = static_cast<size_t>(dims[2]);
  }

  if (rank == 0) {
    RunMaster(size);
  } else {
    RunWorker();
  }

  auto &out = GetOutput();
  if (rank != 0) {
    out.assign(a_rows_ * b_cols_, 0.0);
  }

  MPI_Bcast(out.data(), static_cast<int>(out.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return true;
}

void TabalaevAMatrixMulStrassenALL::RunMaster(int size) {
  if (size == 1 || padded_n_ <= kBaseCaseSize) {
    MasterBase(size);
  } else {
    MasterAll(size);
  }

  auto &out = GetOutput();
  out.assign(a_rows_ * b_cols_, 0.0);

  const auto &result_c = result_c_;
  size_t a_rows = a_rows_;
  size_t b_cols = b_cols_;
  size_t padded_n = padded_n_;

#pragma omp parallel for default(none) shared(out, result_c, a_rows, b_cols, padded_n)
  for (size_t i = 0; i < a_rows; ++i) {
    for (size_t j = 0; j < b_cols; ++j) {
      out[(i * b_cols) + j] = result_c[(i * padded_n) + j];
    }
  }
}

void TabalaevAMatrixMulStrassenALL::MasterBase(int size) {
  result_c_ = StrassenMultiply(padded_a_, padded_b_, padded_n_);

  uint64_t terminate_signal = 0;
  for (int dest = 1; dest < size; ++dest) {
    MPI_Send(&terminate_signal, 1, MPI_UINT64_T, dest, 0, MPI_COMM_WORLD);
  }
}

void TabalaevAMatrixMulStrassenALL::MasterAll(int size) {
  size_t h = padded_n_ / 2;

  std::vector<double> a11;
  std::vector<double> a12;
  std::vector<double> a21;
  std::vector<double> a22;

  std::vector<double> b11;
  std::vector<double> b12;
  std::vector<double> b21;
  std::vector<double> b22;

  SplitMatrix(padded_a_, padded_n_, a11, a12, a21, a22);
  SplitMatrix(padded_b_, padded_n_, b11, b12, b21, b22);

  std::vector<std::vector<double>> task_a = {Add(a11, a22),      Add(a21, a22),     a11, a22, Add(a11, a12),
                                             Subtract(a21, a11), Subtract(a12, a22)};
  std::vector<std::vector<double>> task_b = {Add(b11, b22), b11,           Subtract(b12, b22), Subtract(b21, b11),
                                             b22,           Add(b11, b12), Add(b21, b22)};

  std::vector<std::vector<double>> p(7, std::vector<double>(h * h));
  int num_tasks = 7;
  int tasks_sent = 0;
  int tasks_completed = 0;

  int h_squared = static_cast<int>(h * h);

  for (int dest = 1; dest < size && tasks_sent < num_tasks; ++dest) {
    auto h_msg = static_cast<uint64_t>(h);
    MPI_Send(&h_msg, 1, MPI_UINT64_T, dest, 0, MPI_COMM_WORLD);
    MPI_Send(&tasks_sent, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
    MPI_Send(task_a[tasks_sent].data(), h_squared, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
    MPI_Send(task_b[tasks_sent].data(), h_squared, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD);
    tasks_sent++;
  }

  while (tasks_completed < num_tasks) {
    MPI_Status status;
    int task_id = 0;
    MPI_Recv(&task_id, 1, MPI_INT, MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, &status);

    int worker = status.MPI_SOURCE;
    MPI_Recv(p[task_id].data(), h_squared, MPI_DOUBLE, worker, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    tasks_completed++;

    if (tasks_sent < num_tasks) {
      auto h_msg = static_cast<uint64_t>(h);
      MPI_Send(&h_msg, 1, MPI_UINT64_T, worker, 0, MPI_COMM_WORLD);
      MPI_Send(&tasks_sent, 1, MPI_INT, worker, 1, MPI_COMM_WORLD);
      MPI_Send(task_a[tasks_sent].data(), h_squared, MPI_DOUBLE, worker, 2, MPI_COMM_WORLD);
      MPI_Send(task_b[tasks_sent].data(), h_squared, MPI_DOUBLE, worker, 3, MPI_COMM_WORLD);
      tasks_sent++;
    }
  }

  uint64_t terminate_signal = 0;
  for (int dest = 1; dest < size; ++dest) {
    MPI_Send(&terminate_signal, 1, MPI_UINT64_T, dest, 0, MPI_COMM_WORLD);
  }

  std::vector<double> c11(h * h);
  std::vector<double> c12(h * h);
  std::vector<double> c21(h * h);
  std::vector<double> c22(h * h);

#pragma omp parallel for default(none) shared(p, c11, c12, c21, c22, h)
  for (size_t i = 0; i < h * h; ++i) {
    c11[i] = p[0][i] + p[3][i] - p[4][i] + p[6][i];
    c12[i] = p[2][i] + p[4][i];
    c21[i] = p[1][i] + p[3][i];
    c22[i] = p[0][i] - p[1][i] + p[2][i] + p[5][i];
  }

  result_c_ = CombineMatrix(c11, c12, c21, c22, padded_n_);
}

void TabalaevAMatrixMulStrassenALL::RunWorker() {
  while (true) {
    uint64_t h_msg = 0;
    MPI_Status status;
    MPI_Recv(&h_msg, 1, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD, &status);

    if (h_msg == 0) {
      break;
    }
    auto h = static_cast<size_t>(h_msg);
    int h_squared = static_cast<int>(h * h);

    int task_id = 0;
    MPI_Recv(&task_id, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::vector<double> a_sub(h * h);
    std::vector<double> b_sub(h * h);
    MPI_Recv(a_sub.data(), h_squared, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(b_sub.data(), h_squared, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::vector<double> res = StrassenMultiply(a_sub, b_sub, h);

    MPI_Send(&task_id, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
    MPI_Send(res.data(), h_squared, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
  }
}

bool TabalaevAMatrixMulStrassenALL::PostProcessingImpl() {
  return true;
}

std::vector<double> TabalaevAMatrixMulStrassenALL::Add(const std::vector<double> &mat_a,
                                                       const std::vector<double> &mat_b) {
  const size_t n = mat_a.size();
  std::vector<double> res(n);

#pragma omp parallel for default(none) shared(mat_a, mat_b, res, n)
  for (size_t i = 0; i < n; ++i) {
    res[i] = mat_a[i] + mat_b[i];
  }

  return res;
}

std::vector<double> TabalaevAMatrixMulStrassenALL::Subtract(const std::vector<double> &mat_a,
                                                            const std::vector<double> &mat_b) {
  const size_t n = mat_a.size();
  std::vector<double> res(n);

#pragma omp parallel for default(none) shared(mat_a, mat_b, res, n)
  for (size_t i = 0; i < n; ++i) {
    res[i] = mat_a[i] - mat_b[i];
  }

  return res;
}

std::vector<double> TabalaevAMatrixMulStrassenALL::BaseMultiply(const std::vector<double> &mat_a,
                                                                const std::vector<double> &mat_b, size_t n) {
  std::vector<double> res(n * n, 0.0);

#pragma omp parallel for default(none) shared(mat_a, mat_b, res, n)
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      double temp = mat_a[(i * n) + k];
      if (temp == 0.0) {
        continue;
      }
      for (size_t j = 0; j < n; ++j) {
        res[(i * n) + j] += temp * mat_b[(k * n) + j];
      }
    }
  }

  return res;
}

void TabalaevAMatrixMulStrassenALL::SplitMatrix(const std::vector<double> &src, size_t n, std::vector<double> &c11,
                                                std::vector<double> &c12, std::vector<double> &c21,
                                                std::vector<double> &c22) {
  size_t h = n / 2;
  size_t sz = h * h;
  c11.resize(sz);
  c12.resize(sz);
  c21.resize(sz);
  c22.resize(sz);

#pragma omp parallel for collapse(2) default(none) shared(src, c11, c12, c21, c22, h, n)
  for (size_t i = 0; i < h; ++i) {
    for (size_t j = 0; j < h; ++j) {
      size_t src_idx = (i * n) + j;
      size_t dst_idx = (i * h) + j;
      c11[dst_idx] = src[src_idx];
      c12[dst_idx] = src[src_idx + h];
      c21[dst_idx] = src[src_idx + (h * n)];
      c22[dst_idx] = src[src_idx + (h * n) + h];
    }
  }
}

std::vector<double> TabalaevAMatrixMulStrassenALL::CombineMatrix(const std::vector<double> &c11,
                                                                 const std::vector<double> &c12,
                                                                 const std::vector<double> &c21,
                                                                 const std::vector<double> &c22, size_t n) {
  size_t h = n / 2;
  std::vector<double> res(n * n);

#pragma omp parallel for collapse(2) default(none) shared(res, c11, c12, c21, c22, h, n)
  for (size_t i = 0; i < h; ++i) {
    for (size_t j = 0; j < h; ++j) {
      size_t src_idx = (i * h) + j;
      res[(i * n) + j] = c11[src_idx];
      res[(i * n) + j + h] = c12[src_idx];
      res[((i + h) * n) + j] = c21[src_idx];
      res[((i + h) * n) + j + h] = c22[src_idx];
    }
  }
  return res;
}

std::vector<double> TabalaevAMatrixMulStrassenALL::StrassenMultiply(const std::vector<double> &mat_a,
                                                                    const std::vector<double> &mat_b, size_t n) {
  std::stack<StrassenFrameALL> frames;
  std::stack<std::vector<double>> results;

  frames.push({mat_a, mat_b, n, 0});

  while (!frames.empty()) {
    StrassenFrameALL current = std::move(frames.top());
    frames.pop();

    if (current.n <= kBaseCaseSize) {
      results.push(BaseMultiply(current.mat_a, current.mat_b, current.n));
      continue;
    }

    if (current.stage == 8) {
      std::vector<std::vector<double>> p(7);

      for (int i = 6; i >= 0; --i) {
        p[i] = std::move(results.top());
        results.pop();
      }

      size_t h = current.n / 2;
      size_t sz = h * h;
      std::vector<double> c11(sz);
      std::vector<double> c12(sz);
      std::vector<double> c21(sz);
      std::vector<double> c22(sz);

#pragma omp parallel for default(none) shared(p, c11, c12, c21, c22, sz)
      for (size_t i = 0; i < sz; ++i) {
        c11[i] = p[0][i] + p[3][i] - p[4][i] + p[6][i];
        c12[i] = p[2][i] + p[4][i];
        c21[i] = p[1][i] + p[3][i];
        c22[i] = p[0][i] - p[1][i] + p[2][i] + p[5][i];
      }

      results.push(CombineMatrix(c11, c12, c21, c22, current.n));
    } else {
      size_t h = current.n / 2;
      std::vector<double> a11;
      std::vector<double> a12;
      std::vector<double> a21;
      std::vector<double> a22;
      std::vector<double> b11;
      std::vector<double> b12;
      std::vector<double> b21;
      std::vector<double> b22;

      SplitMatrix(current.mat_a, current.n, a11, a12, a21, a22);
      SplitMatrix(current.mat_b, current.n, b11, b12, b21, b22);

      frames.push({{}, {}, current.n, 8});

      frames.push({Subtract(a12, a22), Add(b21, b22), h, 0});
      frames.push({Subtract(a21, a11), Add(b11, b12), h, 0});
      frames.push({Add(a11, a12), b22, h, 0});
      frames.push({a22, Subtract(b21, b11), h, 0});
      frames.push({a11, Subtract(b12, b22), h, 0});
      frames.push({Add(a21, a22), b11, h, 0});
      frames.push({Add(a11, a22), Add(b11, b22), h, 0});
    }
  }

  return std::move(results.top());
}

}  // namespace tabalaev_a_matrix_mul_strassen
