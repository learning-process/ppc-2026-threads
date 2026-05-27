#include "zorin_d_strassen_alg_matrix_seq/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "util/include/util.hpp"
#include "zorin_d_strassen_alg_matrix_seq/common/include/common.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

constexpr std::size_t kCutoff = 128;
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

    const double *m1_row = m1.data() + (i * half);
    const double *m2_row = m2.data() + (i * half);
    const double *m3_row = m3.data() + (i * half);
    const double *m4_row = m4.data() + (i * half);
    const double *m5_row = m5.data() + (i * half);
    const double *m6_row = m6.data() + (i * half);
    const double *m7_row = m7.data() + (i * half);

    for (std::size_t j = 0; j < half; ++j) {
      c11[j] = m1_row[j] + m4_row[j] - m5_row[j] + m7_row[j];
      c12[j] = m3_row[j] + m5_row[j];
      c21[j] = m2_row[j] + m4_row[j];
      c22[j] = m1_row[j] - m2_row[j] + m3_row[j] + m6_row[j];
    }
  }
}

using StrassenFn = void (*)(const double *, std::size_t, const double *, std::size_t, double *, std::size_t,
                            std::size_t);

void StrassenSeqImpl(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                     std::size_t c_stride, std::size_t n);

constexpr StrassenFn kStrassenSeqFn = &StrassenSeqImpl;

void ComputeProduct(const double *a1, std::size_t a1_stride, const double *a2, std::size_t a2_stride, double a2_coeff,
                    const double *b1, std::size_t b1_stride, const double *b2, std::size_t b2_stride, double b2_coeff,
                    std::vector<double> &out, std::size_t n) {
  std::vector<double> lhs(n * n);
  std::vector<double> rhs(n * n);
  AddToBuffer(a1, a1_stride, a2, a2_stride, lhs.data(), n, a2_coeff);
  AddToBuffer(b1, b1_stride, b2, b2_stride, rhs.data(), n, b2_coeff);
  out.assign(n * n, 0.0);
  kStrassenSeqFn(lhs.data(), n, rhs.data(), n, out.data(), n, n);
}

void ComputeProductSingle(const double *a, std::size_t a_stride, const double *b1, std::size_t b1_stride,
                          const double *b2, std::size_t b2_stride, double b2_coeff, std::vector<double> &out,
                          std::size_t n) {
  std::vector<double> rhs(n * n);
  AddToBuffer(b1, b1_stride, b2, b2_stride, rhs.data(), n, b2_coeff);
  out.assign(n * n, 0.0);
  kStrassenSeqFn(a, a_stride, rhs.data(), n, out.data(), n, n);
}

void ComputeProductSingleLeft(const double *a1, std::size_t a1_stride, const double *a2, std::size_t a2_stride,
                              double a2_coeff, const double *b, std::size_t b_stride, std::vector<double> &out,
                              std::size_t n) {
  std::vector<double> lhs(n * n);
  AddToBuffer(a1, a1_stride, a2, a2_stride, lhs.data(), n, a2_coeff);
  out.assign(n * n, 0.0);
  kStrassenSeqFn(lhs.data(), n, b, b_stride, out.data(), n, n);
}

void ComputeProductOmp(const double *a1, std::size_t a1_stride, const double *a2, std::size_t a2_stride,
                       double a2_coeff, const double *b1, std::size_t b1_stride, const double *b2,
                       std::size_t b2_stride, double b2_coeff, std::vector<double> &out, std::size_t n);

void ComputeProductSingleOmp(const double *a, std::size_t a_stride, const double *b1, std::size_t b1_stride,
                             const double *b2, std::size_t b2_stride, double b2_coeff, std::vector<double> &out,
                             std::size_t n);

void ComputeProductSingleLeftOmp(const double *a1, std::size_t a1_stride, const double *a2, std::size_t a2_stride,
                                 double a2_coeff, const double *b, std::size_t b_stride, std::vector<double> &out,
                                 std::size_t n);

void StrassenOmpLocal(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                      std::size_t c_stride, std::size_t n) {
  if (n <= kCutoff || ppc::util::GetNumThreads() <= 1) {
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
      ComputeProduct(a11, a_stride, a22, a_stride, 1.0, b11, b_stride, b22, b_stride, 1.0, m1, half);
#pragma omp task default(none) shared(m2, a21, a22, b11, a_stride, b_stride, half)
      ComputeProductSingleLeft(a21, a_stride, a22, a_stride, 1.0, b11, b_stride, m2, half);
#pragma omp task default(none) shared(m3, a11, b12, b22, a_stride, b_stride, half)
      ComputeProductSingle(a11, a_stride, b12, b_stride, b22, b_stride, -1.0, m3, half);
#pragma omp task default(none) shared(m4, a22, b21, b11, a_stride, b_stride, half)
      ComputeProductSingle(a22, a_stride, b21, b_stride, b11, b_stride, -1.0, m4, half);
#pragma omp task default(none) shared(m5, a11, a12, b22, a_stride, b_stride, half)
      ComputeProductSingleLeft(a11, a_stride, a12, a_stride, 1.0, b22, b_stride, m5, half);
#pragma omp task default(none) shared(m6, a21, a11, b11, b12, a_stride, b_stride, half)
      ComputeProduct(a21, a_stride, a11, a_stride, -1.0, b11, b_stride, b12, b_stride, 1.0, m6, half);
#pragma omp task default(none) shared(m7, a12, a22, b21, b22, a_stride, b_stride, half)
      ComputeProduct(a12, a_stride, a22, a_stride, -1.0, b21, b_stride, b22, b_stride, 1.0, m7, half);
#pragma omp taskwait
    }
  }

  CombineQuadrants(m1, m2, m3, m4, m5, m6, m7, c, c_stride, half);
}

void ComputeProductOmp(const double *a1, std::size_t a1_stride, const double *a2, std::size_t a2_stride,
                       double a2_coeff, const double *b1, std::size_t b1_stride, const double *b2,
                       std::size_t b2_stride, double b2_coeff, std::vector<double> &out, std::size_t n) {
  std::vector<double> lhs(n * n);
  std::vector<double> rhs(n * n);
  AddToBuffer(a1, a1_stride, a2, a2_stride, lhs.data(), n, a2_coeff);
  AddToBuffer(b1, b1_stride, b2, b2_stride, rhs.data(), n, b2_coeff);
  out.assign(n * n, 0.0);
  StrassenOmpLocal(lhs.data(), n, rhs.data(), n, out.data(), n, n);
}

void ComputeProductSingleOmp(const double *a, std::size_t a_stride, const double *b1, std::size_t b1_stride,
                             const double *b2, std::size_t b2_stride, double b2_coeff, std::vector<double> &out,
                             std::size_t n) {
  std::vector<double> rhs(n * n);
  AddToBuffer(b1, b1_stride, b2, b2_stride, rhs.data(), n, b2_coeff);
  out.assign(n * n, 0.0);
  StrassenOmpLocal(a, a_stride, rhs.data(), n, out.data(), n, n);
}

void ComputeProductSingleLeftOmp(const double *a1, std::size_t a1_stride, const double *a2, std::size_t a2_stride,
                                 double a2_coeff, const double *b, std::size_t b_stride, std::vector<double> &out,
                                 std::size_t n) {
  std::vector<double> lhs(n * n);
  AddToBuffer(a1, a1_stride, a2, a2_stride, lhs.data(), n, a2_coeff);
  out.assign(n * n, 0.0);
  StrassenOmpLocal(lhs.data(), n, b, b_stride, out.data(), n, n);
}

void StrassenSeqImpl(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
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

  std::vector<double> m1(half * half);
  std::vector<double> m2(half * half);
  std::vector<double> m3(half * half);
  std::vector<double> m4(half * half);
  std::vector<double> m5(half * half);
  std::vector<double> m6(half * half);
  std::vector<double> m7(half * half);
  std::vector<double> lhs(half * half);
  std::vector<double> rhs(half * half);

  AddToBuffer(a11, a_stride, a22, a_stride, lhs.data(), half, 1.0);
  AddToBuffer(b11, b_stride, b22, b_stride, rhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, rhs.data(), half, m1.data(), half, half);

  AddToBuffer(a21, a_stride, a22, a_stride, lhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, b11, b_stride, m2.data(), half, half);

  AddToBuffer(b12, b_stride, b22, b_stride, rhs.data(), half, -1.0);
  kStrassenSeqFn(a11, a_stride, rhs.data(), half, m3.data(), half, half);

  AddToBuffer(b21, b_stride, b11, b_stride, rhs.data(), half, -1.0);
  kStrassenSeqFn(a22, a_stride, rhs.data(), half, m4.data(), half, half);

  AddToBuffer(a11, a_stride, a12, a_stride, lhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, b22, b_stride, m5.data(), half, half);

  AddToBuffer(a21, a_stride, a11, a_stride, lhs.data(), half, -1.0);
  AddToBuffer(b11, b_stride, b12, b_stride, rhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, rhs.data(), half, m6.data(), half, half);

  AddToBuffer(a12, a_stride, a22, a_stride, lhs.data(), half, -1.0);
  AddToBuffer(b21, b_stride, b22, b_stride, rhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, rhs.data(), half, m7.data(), half, half);

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

  switch (task_id) {
    case 0:
      ComputeProductOmp(a11, a_stride, a22, a_stride, 1.0, b11, b_stride, b22, b_stride, 1.0, m, half);
      AddContribution(local_accum.data(), half * 2, m, 0, 0, half, 1.0);
      AddContribution(local_accum.data(), half * 2, m, half, half, half, 1.0);
      break;
    case 1:
      ComputeProductSingleLeftOmp(a21, a_stride, a22, a_stride, 1.0, b11, b_stride, m, half);
      AddContribution(local_accum.data(), half * 2, m, half, 0, half, 1.0);
      AddContribution(local_accum.data(), half * 2, m, half, half, half, -1.0);
      break;
    case 2:
      ComputeProductSingleOmp(a11, a_stride, b12, b_stride, b22, b_stride, -1.0, m, half);
      AddContribution(local_accum.data(), half * 2, m, 0, half, half, 1.0);
      AddContribution(local_accum.data(), half * 2, m, half, half, half, 1.0);
      break;
    case 3:
      ComputeProductSingleOmp(a22, a_stride, b21, b_stride, b11, b_stride, -1.0, m, half);
      AddContribution(local_accum.data(), half * 2, m, 0, 0, half, 1.0);
      AddContribution(local_accum.data(), half * 2, m, half, 0, half, 1.0);
      break;
    case 4:
      ComputeProductSingleLeftOmp(a11, a_stride, a12, a_stride, 1.0, b22, b_stride, m, half);
      AddContribution(local_accum.data(), half * 2, m, 0, 0, half, -1.0);
      AddContribution(local_accum.data(), half * 2, m, 0, half, half, 1.0);
      break;
    case 5:
      ComputeProductOmp(a21, a_stride, a11, a_stride, -1.0, b11, b_stride, b12, b_stride, 1.0, m, half);
      AddContribution(local_accum.data(), half * 2, m, half, half, half, 1.0);
      break;
    case 6:
      ComputeProductOmp(a12, a_stride, a22, a_stride, -1.0, b21, b_stride, b22, b_stride, 1.0, m, half);
      AddContribution(local_accum.data(), half * 2, m, 0, 0, half, 1.0);
      break;
    default:
      break;
  }
}

}  // namespace

ZorinDStrassenAlgMatrixALL::ZorinDStrassenAlgMatrixALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ZorinDStrassenAlgMatrixALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }

  const auto &in = GetInput();
  if (in.n == 0) {
    return false;
  }
  if (in.a.size() != in.n * in.n) {
    return false;
  }
  if (in.b.size() != in.n * in.n) {
    return false;
  }
  if (!GetOutput().empty()) {
    return false;
  }
  return true;
}

bool ZorinDStrassenAlgMatrixALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    const auto n = GetInput().n;
    GetOutput().assign(n * n, 0.0);
  } else {
    GetOutput().clear();
  }
  return true;
}

bool ZorinDStrassenAlgMatrixALL::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  std::uint64_t n_u64 = 0;
  if (rank == 0) {
    n_u64 = static_cast<std::uint64_t>(GetInput().n);
  }
  MPI_Bcast(&n_u64, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  const auto n = static_cast<std::size_t>(n_u64);
  const std::size_t padded = NextPow2(n);

  std::vector<double> a_pad(padded * padded, 0.0);
  std::vector<double> b_pad(padded * padded, 0.0);

  if (rank == 0) {
    const auto &in = GetInput();
    for (std::size_t i = 0; i < n; ++i) {
      std::copy_n(in.a.data() + (i * n), n, a_pad.data() + (i * padded));
      std::copy_n(in.b.data() + (i * n), n, b_pad.data() + (i * padded));
    }
  }

  MPI_Bcast(a_pad.data(), static_cast<int>(a_pad.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(b_pad.data(), static_cast<int>(b_pad.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::vector<double> global_c(padded * padded, 0.0);

  if (padded <= kCutoff || world_size == 1) {
    if (rank == 0) {
      StrassenOmpLocal(a_pad.data(), padded, b_pad.data(), padded, global_c.data(), padded, padded);
    }
  } else {
    const std::size_t half = padded / 2;
    std::vector<double> local_c(padded * padded, 0.0);

    for (int task_id = rank; task_id < 7; task_id += world_size) {
      ComputeAssignedProduct(task_id, a_pad.data(), padded, b_pad.data(), padded, half, local_c);
    }

    MPI_Reduce(local_c.data(), global_c.data(), static_cast<int>(global_c.size()), MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
  }

  auto &out = GetOutput();
  if (rank != 0) {
    out.assign(n * n, 0.0);
  }
  if (rank == 0) {
    for (std::size_t i = 0; i < n; ++i) {
      std::copy_n(global_c.data() + (i * padded), n, out.data() + (i * n));
    }
  }

  MPI_Bcast(out.data(), static_cast<int>(out.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return true;
}

bool ZorinDStrassenAlgMatrixALL::PostProcessingImpl() {
  return true;
}

}  // namespace zorin_d_strassen_alg_matrix_seq
