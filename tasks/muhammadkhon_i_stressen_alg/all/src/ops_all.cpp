#include "muhammadkhon_i_stressen_alg/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include "muhammadkhon_i_stressen_alg/common/include/common.hpp"
#include "util/include/util.hpp"

namespace muhammadkhon_i_stressen_alg {

namespace {

constexpr std::size_t kCutoff = 64;
constexpr std::size_t kBlockSize = 64;

std::size_t NextPow2(std::size_t x) {
  if (x <= 1) {
    return 1;
  }
  std::size_t p = 1;
  while (p < x) {
    p <<= 1;
  }
  return p;
}

void ZeroMatrix(double *dst, std::size_t stride, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    std::fill_n(dst + (i * stride), n, 0.0);
  }
}

void AddToBuffer(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *dst,
                 std::size_t n, double b_coeff) {
  for (std::size_t i = 0; i < n; ++i) {
    const double *a_row = a + (i * a_stride);
    const double *b_row = b + (i * b_stride);
    double *dst_row = dst + (i * n);
    for (std::size_t j = 0; j < n; ++j) {
      dst_row[j] = a_row[j] + (b_coeff * b_row[j]);
    }
  }
}

void MulMicroBlock(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                   std::size_t c_stride, std::size_t i_begin, std::size_t i_end, std::size_t k_begin, std::size_t k_end,
                   std::size_t j_begin, std::size_t j_end) {
  for (std::size_t i = i_begin; i < i_end; ++i) {
    double *c_row = c + (i * c_stride);
    const double *a_row = a + (i * a_stride);
    for (std::size_t k = k_begin; k < k_end; ++k) {
      const double aik = a_row[k];
      const double *b_row = b + (k * b_stride);
      for (std::size_t j = j_begin; j < j_end; ++j) {
        c_row[j] += aik * b_row[j];
      }
    }
  }
}

void NaiveMulBlocked(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                     std::size_t c_stride, std::size_t n) {
  ZeroMatrix(c, c_stride, n);

  const auto n_signed = static_cast<std::ptrdiff_t>(n);
  const auto block_signed = static_cast<std::ptrdiff_t>(kBlockSize);

#pragma omp parallel for schedule(static) default(none) \
    shared(a, a_stride, b, b_stride, c, c_stride, n, n_signed, block_signed)
  for (std::ptrdiff_t ii = 0; ii < n_signed; ii += block_signed) {
    const auto ii_usize = static_cast<std::size_t>(ii);
    const std::size_t i_end = std::min(ii_usize + kBlockSize, n);
    for (std::size_t kk = 0; kk < n; kk += kBlockSize) {
      const std::size_t k_end = std::min(kk + kBlockSize, n);
      for (std::size_t jj = 0; jj < n; jj += kBlockSize) {
        const std::size_t j_end = std::min(jj + kBlockSize, n);
        MulMicroBlock(a, a_stride, b, b_stride, c, c_stride, ii_usize, i_end, kk, k_end, jj, j_end);
      }
    }
  }
}

void CombineQuadrants(const std::vector<double> &m1, const std::vector<double> &m2, const std::vector<double> &m3,
                      const std::vector<double> &m4, const std::vector<double> &m5, const std::vector<double> &m6,
                      const std::vector<double> &m7, double *c, std::size_t c_stride, std::size_t half) {
  for (std::size_t i = 0; i < half; ++i) {
    double *c11 = c + (i * c_stride);
    double *c12 = c11 + half;
    double *c21 = c + ((i + half) * c_stride);
    double *c22 = c21 + half;
    const double *m1r = m1.data() + (i * half);
    const double *m2r = m2.data() + (i * half);
    const double *m3r = m3.data() + (i * half);
    const double *m4r = m4.data() + (i * half);
    const double *m5r = m5.data() + (i * half);
    const double *m6r = m6.data() + (i * half);
    const double *m7r = m7.data() + (i * half);
    for (std::size_t j = 0; j < half; ++j) {
      c11[j] = m1r[j] + m4r[j] - m5r[j] + m7r[j];
      c12[j] = m3r[j] + m5r[j];
      c21[j] = m2r[j] + m4r[j];
      c22[j] = m1r[j] - m2r[j] + m3r[j] + m6r[j];
    }
  }
}

// Последовательный Штрассен через std::function (без рекурсивных свободных функций)
void StrassenSeq(const double *a_in, std::size_t a_stride_in, const double *b_in, std::size_t b_stride_in, double *c_in,
                 std::size_t c_stride_in, std::size_t n_in) {
  std::function<void(const double *, std::size_t, const double *, std::size_t, double *, std::size_t, std::size_t)>
      impl = [&](const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                 std::size_t c_stride, std::size_t n) {
    if (n <= kCutoff) {
      NaiveMulBlocked(a, a_stride, b, b_stride, c, c_stride, n);
      return;
    }
    const std::size_t half = n / 2;

    const double *a11 = a;
    const double *a12 = a + half;
    const double *a21 = a + (half * a_stride);
    const double *a22 = a21 + half;
    const double *b11 = b;
    const double *b12 = b + half;
    const double *b21 = b + (half * b_stride);
    const double *b22 = b21 + half;

    std::vector<double> lhs(half * half);
    std::vector<double> rhs(half * half);
    std::vector<double> m1(half * half);
    std::vector<double> m2(half * half);
    std::vector<double> m3(half * half);
    std::vector<double> m4(half * half);
    std::vector<double> m5(half * half);
    std::vector<double> m6(half * half);
    std::vector<double> m7(half * half);

    AddToBuffer(a11, a_stride, a22, a_stride, lhs.data(), half, 1.0);
    AddToBuffer(b11, b_stride, b22, b_stride, rhs.data(), half, 1.0);
    impl(lhs.data(), half, rhs.data(), half, m1.data(), half, half);

    AddToBuffer(a21, a_stride, a22, a_stride, lhs.data(), half, 1.0);
    impl(lhs.data(), half, b11, b_stride, m2.data(), half, half);

    AddToBuffer(b12, b_stride, b22, b_stride, rhs.data(), half, -1.0);
    impl(a11, a_stride, rhs.data(), half, m3.data(), half, half);

    AddToBuffer(b21, b_stride, b11, b_stride, rhs.data(), half, -1.0);
    impl(a22, a_stride, rhs.data(), half, m4.data(), half, half);

    AddToBuffer(a11, a_stride, a12, a_stride, lhs.data(), half, 1.0);
    impl(lhs.data(), half, b22, b_stride, m5.data(), half, half);

    AddToBuffer(a21, a_stride, a11, a_stride, lhs.data(), half, -1.0);
    AddToBuffer(b11, b_stride, b12, b_stride, rhs.data(), half, 1.0);
    impl(lhs.data(), half, rhs.data(), half, m6.data(), half, half);

    AddToBuffer(a12, a_stride, a22, a_stride, lhs.data(), half, -1.0);
    AddToBuffer(b21, b_stride, b22, b_stride, rhs.data(), half, 1.0);
    impl(lhs.data(), half, rhs.data(), half, m7.data(), half, half);

    CombineQuadrants(m1, m2, m3, m4, m5, m6, m7, c, c_stride, half);
  };

  impl(a_in, a_stride_in, b_in, b_stride_in, c_in, c_stride_in, n_in);
}

// OMP-параллельный Штрассен (верхний уровень через tasks, базовый — OMP parallel for)
void StrassenOmpLocal(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                      std::size_t c_stride, std::size_t n) {
  if (n <= kCutoff || ppc::util::GetNumThreads() <= 1) {
    StrassenSeq(a, a_stride, b, b_stride, c, c_stride, n);
    return;
  }

  const std::size_t half = n / 2;

  const double *a11 = a;
  const double *a12 = a + half;
  const double *a21 = a + (half * a_stride);
  const double *a22 = a21 + half;
  const double *b11 = b;
  const double *b12 = b + half;
  const double *b21 = b + (half * b_stride);
  const double *b22 = b21 + half;

  std::vector<double> m1;
  std::vector<double> m2;
  std::vector<double> m3;
  std::vector<double> m4;
  std::vector<double> m5;
  std::vector<double> m6;
  std::vector<double> m7;

#pragma omp parallel default(none) \
    shared(m1, m2, m3, m4, m5, m6, m7, a11, a12, a21, a22, b11, b12, b21, b22, a_stride, b_stride, half)
  {
#pragma omp single nowait
    {
#pragma omp task default(none) shared(m1, a11, a22, b11, b22, a_stride, b_stride, half)
      {
        std::vector<double> lhs(half * half);
        std::vector<double> rhs(half * half);
        AddToBuffer(a11, a_stride, a22, a_stride, lhs.data(), half, 1.0);
        AddToBuffer(b11, b_stride, b22, b_stride, rhs.data(), half, 1.0);
        m1.assign(half * half, 0.0);
        StrassenSeq(lhs.data(), half, rhs.data(), half, m1.data(), half, half);
      }
#pragma omp task default(none) shared(m2, a21, a22, b11, a_stride, b_stride, half)
      {
        std::vector<double> lhs(half * half);
        AddToBuffer(a21, a_stride, a22, a_stride, lhs.data(), half, 1.0);
        m2.assign(half * half, 0.0);
        StrassenSeq(lhs.data(), half, b11, b_stride, m2.data(), half, half);
      }
#pragma omp task default(none) shared(m3, a11, b12, b22, a_stride, b_stride, half)
      {
        std::vector<double> rhs(half * half);
        AddToBuffer(b12, b_stride, b22, b_stride, rhs.data(), half, -1.0);
        m3.assign(half * half, 0.0);
        StrassenSeq(a11, a_stride, rhs.data(), half, m3.data(), half, half);
      }
#pragma omp task default(none) shared(m4, a22, b21, b11, a_stride, b_stride, half)
      {
        std::vector<double> rhs(half * half);
        AddToBuffer(b21, b_stride, b11, b_stride, rhs.data(), half, -1.0);
        m4.assign(half * half, 0.0);
        StrassenSeq(a22, a_stride, rhs.data(), half, m4.data(), half, half);
      }
#pragma omp task default(none) shared(m5, a11, a12, b22, a_stride, b_stride, half)
      {
        std::vector<double> lhs(half * half);
        AddToBuffer(a11, a_stride, a12, a_stride, lhs.data(), half, 1.0);
        m5.assign(half * half, 0.0);
        StrassenSeq(lhs.data(), half, b22, b_stride, m5.data(), half, half);
      }
#pragma omp task default(none) shared(m6, a21, a11, b11, b12, a_stride, b_stride, half)
      {
        std::vector<double> lhs(half * half);
        std::vector<double> rhs(half * half);
        AddToBuffer(a21, a_stride, a11, a_stride, lhs.data(), half, -1.0);
        AddToBuffer(b11, b_stride, b12, b_stride, rhs.data(), half, 1.0);
        m6.assign(half * half, 0.0);
        StrassenSeq(lhs.data(), half, rhs.data(), half, m6.data(), half, half);
      }
#pragma omp task default(none) shared(m7, a12, a22, b21, b22, a_stride, b_stride, half)
      {
        std::vector<double> lhs(half * half);
        std::vector<double> rhs(half * half);
        AddToBuffer(a12, a_stride, a22, a_stride, lhs.data(), half, -1.0);
        AddToBuffer(b21, b_stride, b22, b_stride, rhs.data(), half, 1.0);
        m7.assign(half * half, 0.0);
        StrassenSeq(lhs.data(), half, rhs.data(), half, m7.data(), half, half);
      }
#pragma omp taskwait
    }
  }

  CombineQuadrants(m1, m2, m3, m4, m5, m6, m7, c, c_stride, half);
}

void AddContribution(double *accum, std::size_t stride, const std::vector<double> &block, std::size_t row_offset,
                     std::size_t col_offset, std::size_t half, double coeff) {
  for (std::size_t i = 0; i < half; ++i) {
    double *dst_row = accum + ((row_offset + i) * stride) + col_offset;
    const double *src_row = block.data() + (i * half);
    for (std::size_t j = 0; j < half; ++j) {
      dst_row[j] += coeff * src_row[j];
    }
  }
}

void ComputeAssignedProduct(int task_id, const double *a, std::size_t a_stride, const double *b, std::size_t b_stride,
                            std::size_t half, std::vector<double> &local_accum) {
  const double *a11 = a;
  const double *a12 = a + half;
  const double *a21 = a + (half * a_stride);
  const double *a22 = a21 + half;
  const double *b11 = b;
  const double *b12 = b + half;
  const double *b21 = b + (half * b_stride);
  const double *b22 = b21 + half;

  std::vector<double> m(half * half, 0.0);
  std::vector<double> lhs(half * half);
  std::vector<double> rhs(half * half);

  switch (task_id) {
    case 0:  // M1 = (A11+A22)(B11+B22) -> C11, C22
      AddToBuffer(a11, a_stride, a22, a_stride, lhs.data(), half, 1.0);
      AddToBuffer(b11, b_stride, b22, b_stride, rhs.data(), half, 1.0);
      StrassenOmpLocal(lhs.data(), half, rhs.data(), half, m.data(), half, half);
      AddContribution(local_accum.data(), half * 2, m, 0, 0, half, 1.0);
      AddContribution(local_accum.data(), half * 2, m, half, half, half, 1.0);
      break;
    case 1:  // M2 = (A21+A22)B11 -> C21, C22
      AddToBuffer(a21, a_stride, a22, a_stride, lhs.data(), half, 1.0);
      StrassenOmpLocal(lhs.data(), half, b11, b_stride, m.data(), half, half);
      AddContribution(local_accum.data(), half * 2, m, half, 0, half, 1.0);
      AddContribution(local_accum.data(), half * 2, m, half, half, half, -1.0);
      break;
    case 2:  // M3 = A11(B12-B22) -> C12, C22
      AddToBuffer(b12, b_stride, b22, b_stride, rhs.data(), half, -1.0);
      StrassenOmpLocal(a11, a_stride, rhs.data(), half, m.data(), half, half);
      AddContribution(local_accum.data(), half * 2, m, 0, half, half, 1.0);
      AddContribution(local_accum.data(), half * 2, m, half, half, half, 1.0);
      break;
    case 3:  // M4 = A22(B21-B11) -> C11, C21
      AddToBuffer(b21, b_stride, b11, b_stride, rhs.data(), half, -1.0);
      StrassenOmpLocal(a22, a_stride, rhs.data(), half, m.data(), half, half);
      AddContribution(local_accum.data(), half * 2, m, 0, 0, half, 1.0);
      AddContribution(local_accum.data(), half * 2, m, half, 0, half, 1.0);
      break;
    case 4:  // M5 = (A11+A12)B22 -> C11, C12
      AddToBuffer(a11, a_stride, a12, a_stride, lhs.data(), half, 1.0);
      StrassenOmpLocal(lhs.data(), half, b22, b_stride, m.data(), half, half);
      AddContribution(local_accum.data(), half * 2, m, 0, 0, half, -1.0);
      AddContribution(local_accum.data(), half * 2, m, 0, half, half, 1.0);
      break;
    case 5:  // M6 = (A21-A11)(B11+B12) -> C22
      AddToBuffer(a21, a_stride, a11, a_stride, lhs.data(), half, -1.0);
      AddToBuffer(b11, b_stride, b12, b_stride, rhs.data(), half, 1.0);
      StrassenOmpLocal(lhs.data(), half, rhs.data(), half, m.data(), half, half);
      AddContribution(local_accum.data(), half * 2, m, half, half, half, 1.0);
      break;
    case 6:  // M7 = (A12-A22)(B21+B22) -> C11
      AddToBuffer(a12, a_stride, a22, a_stride, lhs.data(), half, -1.0);
      AddToBuffer(b21, b_stride, b22, b_stride, rhs.data(), half, 1.0);
      StrassenOmpLocal(lhs.data(), half, rhs.data(), half, m.data(), half, half);
      AddContribution(local_accum.data(), half * 2, m, 0, 0, half, 1.0);
      break;
    default:
      break;
  }
}

}  // namespace

MuhammadkhonIStressenAlgALL::MuhammadkhonIStressenAlgALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MuhammadkhonIStressenAlgALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }
  const auto &in = GetInput();
  return in.a_rows > 0 && in.a_cols_b_rows > 0 && in.b_cols > 0 &&
         in.a.size() == static_cast<size_t>(in.a_rows * in.a_cols_b_rows) &&
         in.b.size() == static_cast<size_t>(in.a_cols_b_rows * in.b_cols);
}

bool MuhammadkhonIStressenAlgALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetOutput() = {};
    const auto &in = GetInput();
    a_rows_ = in.a_rows;
    a_cols_b_rows_ = in.a_cols_b_rows;
    b_cols_ = in.b_cols;

    const size_t max_dim = std::max({a_rows_, a_cols_b_rows_, b_cols_});
    padded_n_ = NextPow2(max_dim);

    padded_a_.assign(padded_n_ * padded_n_, 0.0);
    padded_b_.assign(padded_n_ * padded_n_, 0.0);

    for (size_t i = 0; i < a_rows_; ++i) {
      for (size_t j = 0; j < a_cols_b_rows_; ++j) {
        padded_a_[(i * padded_n_) + j] = in.a[(i * a_cols_b_rows_) + j];
      }
    }
    for (size_t i = 0; i < a_cols_b_rows_; ++i) {
      for (size_t j = 0; j < b_cols_; ++j) {
        padded_b_[(i * padded_n_) + j] = in.b[(i * b_cols_) + j];
      }
    }
  } else {
    GetOutput().clear();
  }
  return true;
}

bool MuhammadkhonIStressenAlgALL::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Broadcast размеры
  std::array<std::uint64_t, 4> dims = {0, 0, 0, 0};
  if (rank == 0) {
    dims[0] = static_cast<std::uint64_t>(a_rows_);
    dims[1] = static_cast<std::uint64_t>(a_cols_b_rows_);
    dims[2] = static_cast<std::uint64_t>(b_cols_);
    dims[3] = static_cast<std::uint64_t>(padded_n_);
  }
  MPI_Bcast(dims.data(), 4, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    a_rows_ = static_cast<size_t>(dims[0]);
    a_cols_b_rows_ = static_cast<size_t>(dims[1]);
    b_cols_ = static_cast<size_t>(dims[2]);
    padded_n_ = static_cast<size_t>(dims[3]);
    padded_a_.assign(padded_n_ * padded_n_, 0.0);
    padded_b_.assign(padded_n_ * padded_n_, 0.0);
  }

  MPI_Bcast(padded_a_.data(), static_cast<int>(padded_a_.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(padded_b_.data(), static_cast<int>(padded_b_.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  result_c_.assign(padded_n_ * padded_n_, 0.0);

  if (padded_n_ <= kCutoff || world_size == 1) {
    if (rank == 0) {
      StrassenOmpLocal(padded_a_.data(), padded_n_, padded_b_.data(), padded_n_, result_c_.data(), padded_n_,
                       padded_n_);
    }
  } else {
    const std::size_t half = padded_n_ / 2;
    std::vector<double> local_c(padded_n_ * padded_n_, 0.0);

    for (int task_id = rank; task_id < 7; task_id += world_size) {
      ComputeAssignedProduct(task_id, padded_a_.data(), padded_n_, padded_b_.data(), padded_n_, half, local_c);
    }

    MPI_Reduce(local_c.data(), result_c_.data(), static_cast<int>(result_c_.size()), MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
  }

  auto &out = GetOutput();
  out.assign(a_rows_ * b_cols_, 0.0);
  if (rank == 0) {
    for (size_t i = 0; i < a_rows_; ++i) {
      for (size_t j = 0; j < b_cols_; ++j) {
        out[(i * b_cols_) + j] = result_c_[(i * padded_n_) + j];
      }
    }
  }
  MPI_Bcast(out.data(), static_cast<int>(out.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return true;
}

bool MuhammadkhonIStressenAlgALL::PostProcessingImpl() {
  return true;
}

}  // namespace muhammadkhon_i_stressen_alg
