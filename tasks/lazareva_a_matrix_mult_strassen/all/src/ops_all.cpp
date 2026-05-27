#include "lazareva_a_matrix_mult_strassen/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "lazareva_a_matrix_mult_strassen/common/include/common.hpp"

namespace lazareva_a_matrix_mult_strassen {

LazarevaATestTaskALL::LazarevaATestTaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LazarevaATestTaskALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int ok = 0;
  if (rank == 0) {
    const auto &input = GetInput();
    if (input.n > 0 && input.a.size() == static_cast<size_t>(input.n) * input.n &&
        input.b.size() == static_cast<size_t>(input.n) * input.n) {
      ok = 1;
    }
  }
  MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return ok == 1;
}

bool LazarevaATestTaskALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    n_ = GetInput().n;
    padded_n_ = NextPowerOfTwo(n_);
    a_ = PadMatrix(GetInput().a, n_, padded_n_);
    b_ = PadMatrix(GetInput().b, n_, padded_n_);
  }
  MPI_Bcast(&n_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&padded_n_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return true;
}

bool LazarevaATestTaskALL::RunImpl() {
  result_ = StrassenALL(a_, b_, padded_n_);
  return true;
}

bool LazarevaATestTaskALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const size_t final_sz = static_cast<size_t>(n_) * n_;

  if (rank == 0) {
    GetOutput() = UnpadMatrix(result_, padded_n_, n_);
  } else {
    GetOutput().assign(final_sz, 0.0);
  }

  MPI_Bcast(GetOutput().data(), static_cast<int>(final_sz), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return true;
}

int LazarevaATestTaskALL::NextPowerOfTwo(int n) {
  int p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

std::vector<double> LazarevaATestTaskALL::PadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  const size_t new_sz = static_cast<size_t>(new_n) * new_n;
  std::vector<double> res(new_sz, 0.0);
  for (int i = 0; i < old_n; ++i) {
    std::copy(m.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * old_n),
              m.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i + 1) * old_n),
              res.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * new_n));
  }
  return res;
}

std::vector<double> LazarevaATestTaskALL::UnpadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  const size_t new_sz = static_cast<size_t>(new_n) * new_n;
  std::vector<double> res(new_sz);
  for (int i = 0; i < new_n; ++i) {
    std::copy(m.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * old_n),
              m.begin() + static_cast<ptrdiff_t>((static_cast<size_t>(i) * old_n) + new_n),
              res.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * new_n));
  }
  return res;
}

std::vector<double> LazarevaATestTaskALL::Add(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const size_t sz = static_cast<size_t>(n) * n;
  std::vector<double> res(sz);
  for (size_t i = 0; i < sz; ++i) {
    res[i] = a[i] + b[i];
  }
  return res;
}

std::vector<double> LazarevaATestTaskALL::Sub(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const size_t sz = static_cast<size_t>(n) * n;
  std::vector<double> res(sz);
  for (size_t i = 0; i < sz; ++i) {
    res[i] = a[i] - b[i];
  }
  return res;
}

void LazarevaATestTaskALL::Split(const std::vector<double> &p, int n, std::vector<double> &a11,
                                 std::vector<double> &a12, std::vector<double> &a21, std::vector<double> &a22) {
  const int h = n / 2;
  const size_t h_sz = static_cast<size_t>(h) * h;
  a11.resize(h_sz);
  a12.resize(h_sz);
  a21.resize(h_sz);
  a22.resize(h_sz);

  for (int i = 0; i < h; ++i) {
    const double *src_top = p.data() + (static_cast<size_t>(i) * n);
    const double *src_bot = p.data() + (static_cast<size_t>(i + h) * n);
    std::copy(src_top, src_top + h, a11.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * h));
    std::copy(src_top + h, src_top + n, a12.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * h));
    std::copy(src_bot, src_bot + h, a21.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * h));
    std::copy(src_bot + h, src_bot + n, a22.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * h));
  }
}

std::vector<double> LazarevaATestTaskALL::Merge(const std::vector<double> &c11, const std::vector<double> &c12,
                                                const std::vector<double> &c21, const std::vector<double> &c22, int h) {
  const int n = h * 2;
  std::vector<double> res(static_cast<size_t>(n) * n);

  for (int i = 0; i < h; ++i) {
    double *dst_top = res.data() + (static_cast<size_t>(i) * n);
    double *dst_bot = res.data() + (static_cast<size_t>(i + h) * n);
    std::copy(c11.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * h),
              c11.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i + 1) * h), dst_top);
    std::copy(c12.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * h),
              c12.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i + 1) * h), dst_top + h);
    std::copy(c21.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * h),
              c21.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i + 1) * h), dst_bot);
    std::copy(c22.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i) * h),
              c22.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(i + 1) * h), dst_bot + h);
  }
  return res;
}

std::vector<double> LazarevaATestTaskALL::NaiveMult(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto n_sz = static_cast<size_t>(n);
  std::vector<double> c(n_sz * n_sz, 0.0);

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const double aik = a[(static_cast<size_t>(i) * n_sz) + static_cast<size_t>(k)];
      if (std::abs(aik) < 1e-18) {
        continue;
      }
      const double *b_row = b.data() + (static_cast<size_t>(k) * n_sz);
      double *c_row = c.data() + (static_cast<size_t>(i) * n_sz);
      for (int j = 0; j < n; ++j) {
        c_row[j] += aik * b_row[j];
      }
    }
  }
  return c;
}

std::vector<double> LazarevaATestTaskALL::StrassenTBB(const std::vector<double> &a, const std::vector<double> &b,
                                                      int n) {
  if (n <= 64) {
    return NaiveMult(a, b, n);
  }

  const int h = n / 2;
  std::vector<double> a11;
  std::vector<double> a12;
  std::vector<double> a21;
  std::vector<double> a22;
  std::vector<double> b11;
  std::vector<double> b12;
  std::vector<double> b21;
  std::vector<double> b22;

  Split(a, n, a11, a12, a21, a22);
  Split(b, n, b11, b12, b21, b22);

  std::vector<double> m0;
  std::vector<double> m1;
  std::vector<double> m2;
  std::vector<double> m3;
  std::vector<double> m4;
  std::vector<double> m5;
  std::vector<double> m6;

  oneapi::tbb::parallel_invoke([&]() { m0 = StrassenTBB(Add(a11, a22, h), Add(b11, b22, h), h); }, [&]() {
    m1 = StrassenTBB(Add(a21, a22, h), b11, h);
  }, [&]() { m2 = StrassenTBB(a11, Sub(b12, b22, h), h); }, [&]() {
    m3 = StrassenTBB(a22, Sub(b21, b11, h), h);
  }, [&]() { m4 = StrassenTBB(Add(a11, a12, h), b22, h); }, [&]() {
    m5 = StrassenTBB(Sub(a21, a11, h), Add(b11, b12, h), h);
  }, [&]() { m6 = StrassenTBB(Sub(a12, a22, h), Add(b21, b22, h), h); });

  const auto c11 = Add(Sub(Add(m0, m3, h), m4, h), m6, h);
  const auto c12 = Add(m2, m4, h);
  const auto c21 = Add(m1, m3, h);
  const auto c22 = Add(Sub(Add(m0, m2, h), m1, h), m5, h);

  return Merge(c11, c12, c21, c22, h);
}

std::vector<double> LazarevaATestTaskALL::StrassenMaster(const std::vector<double> &a, const std::vector<double> &b,
                                                         int h, size_t h_sz, int size) {
  std::vector<double> a11;
  std::vector<double> a12;
  std::vector<double> a21;
  std::vector<double> a22;
  std::vector<double> b11;
  std::vector<double> b12;
  std::vector<double> b21;
  std::vector<double> b22;

  const int n = h * 2;
  Split(a, n, a11, a12, a21, a22);
  Split(b, n, b11, b12, b21, b22);

  std::vector<std::vector<double>> lhs(7);
  std::vector<std::vector<double>> rhs(7);

  lhs[0] = Add(a11, a22, h);
  rhs[0] = Add(b11, b22, h);
  lhs[1] = Add(a21, a22, h);
  rhs[1] = b11;
  lhs[2] = a11;
  rhs[2] = Sub(b12, b22, h);
  lhs[3] = a22;
  rhs[3] = Sub(b21, b11, h);
  lhs[4] = Add(a11, a12, h);
  rhs[4] = b22;
  lhs[5] = Sub(a21, a11, h);
  rhs[5] = Add(b11, b12, h);
  lhs[6] = Sub(a12, a22, h);
  rhs[6] = Add(b21, b22, h);

  std::vector<MPI_Request> send_requests(14, MPI_REQUEST_NULL);
  int req_idx = 0;

  for (int k = 0; k < 7; ++k) {
    const int dest = k % size;
    if (dest != 0) {
      MPI_Isend(lhs[k].data(), static_cast<int>(h_sz), MPI_DOUBLE, dest, k * 2, MPI_COMM_WORLD,
                &send_requests[req_idx]);
      ++req_idx;
      MPI_Isend(rhs[k].data(), static_cast<int>(h_sz), MPI_DOUBLE, dest, (k * 2) + 1, MPI_COMM_WORLD,
                &send_requests[req_idx]);
      ++req_idx;
    }
  }

  std::vector<std::vector<double>> m(7);
  for (int k = 0; k < 7; ++k) {
    if (k % size == 0) {
      m[k] = StrassenTBB(lhs[k], rhs[k], h);
    }
  }

  for (int k = 0; k < 7; ++k) {
    const int src = k % size;
    if (src != 0) {
      m[k].resize(h_sz);
      MPI_Recv(m[k].data(), static_cast<int>(h_sz), MPI_DOUBLE, src, k + 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  if (req_idx > 0) {
    MPI_Waitall(req_idx, send_requests.data(), MPI_STATUSES_IGNORE);
  }

  const auto c11 = Add(Sub(Add(m[0], m[3], h), m[4], h), m[6], h);
  const auto c12 = Add(m[2], m[4], h);
  const auto c21 = Add(m[1], m[3], h);
  const auto c22 = Add(Sub(Add(m[0], m[2], h), m[1], h), m[5], h);

  return Merge(c11, c12, c21, c22, h);
}

void LazarevaATestTaskALL::StrassenWorker(int rank, int h, size_t h_sz, int size) {
  for (int k = 0; k < 7; ++k) {
    if (k % size == rank) {
      std::vector<double> l(h_sz);
      std::vector<double> r(h_sz);
      MPI_Recv(l.data(), static_cast<int>(h_sz), MPI_DOUBLE, 0, k * 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(r.data(), static_cast<int>(h_sz), MPI_DOUBLE, 0, (k * 2) + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      const auto res = StrassenTBB(l, r, h);
      MPI_Send(res.data(), static_cast<int>(h_sz), MPI_DOUBLE, 0, k + 100, MPI_COMM_WORLD);
    }
  }
}

std::vector<double> LazarevaATestTaskALL::StrassenALL(const std::vector<double> &a, const std::vector<double> &b,
                                                      int n) {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (n <= 64 || size == 1) {
    if (rank == 0) {
      return StrassenTBB(a, b, n);
    }
    return {};
  }

  const int h = n / 2;
  const size_t h_sz = static_cast<size_t>(h) * h;

  if (rank == 0) {
    return StrassenMaster(a, b, h, h_sz, size);
  }

  StrassenWorker(rank, h, h_sz, size);
  return {};
}

}  // namespace lazareva_a_matrix_mult_strassen
