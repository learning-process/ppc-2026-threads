#include "smyshlaev_a_sle_cg_seq/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "smyshlaev_a_sle_cg_seq/common/include/common.hpp"
#include "util/include/util.hpp"

namespace smyshlaev_a_sle_cg_seq {

SmyshlaevASleCgTaskALL::SmyshlaevASleCgTaskALL(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  if (ppc::util::IsUnderMpirun()) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  if (rank == 0) {
    GetInput() = in;
  }
}

bool SmyshlaevASleCgTaskALL::ValidationImpl() {
  int rank = 0;
  if (ppc::util::IsUnderMpirun()) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }

  int error = 0;
  if (rank == 0) {
    const auto &a = GetInput().A;
    const auto &b = GetInput().b;
    if (a.empty() || b.empty() || a.size() != b.size() || a.size() != a[0].size()) {
      error = 1;
    }
  }
  if (ppc::util::IsUnderMpirun()) {
    MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  return error == 0;
}

bool SmyshlaevASleCgTaskALL::PreProcessingImpl() {
  int rank = 0;
  if (ppc::util::IsUnderMpirun()) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }

  if (rank == 0) {
    const auto &a = GetInput().A;
    n_ = static_cast<int>(a.size());
    flat_A_.resize(static_cast<size_t>(n_) * n_);
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        flat_A_[(static_cast<size_t>(i) * n_) + j] = a[i][j];
      }
    }
  }
  return true;
}

void SmyshlaevASleCgTaskALL::DistributeInitialData(int rank, bool is_mpi, std::vector<double> &b) {
  if (is_mpi) {
    MPI_Bcast(&n_, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (n_ > 0) {
      if (rank != 0) {
        b.resize(n_);
        flat_A_.resize(static_cast<size_t>(n_) * n_);
      }
      MPI_Bcast(flat_A_.data(), static_cast<int>(flat_A_.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(b.data(), n_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  }
}

double SmyshlaevASleCgTaskALL::ComputeDotProductAll(const std::vector<double> &v1, const std::vector<double> &v2,
                                                    int start, int end, bool is_mpi) {
  double local_sum = 0.0;
#pragma omp parallel for default(none) shared(start, end, v1, v2) reduction(+ : local_sum)
  for (int i = start; i < end; ++i) {
    local_sum += v1[i] * v2[i];
  }
  double global_sum = local_sum;
  if (is_mpi) {
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  return global_sum;
}

void SmyshlaevASleCgTaskALL::ComputeApAll(const std::vector<double> &p, std::vector<double> &ap, int start,
                                          int end) const {
  const auto &m = flat_A_;
  const int l_n = n_;

#pragma omp parallel for default(none) shared(start, end, p, ap, l_n, m)
  for (int i = start; i < end; ++i) {
    double sum = 0.0;
    for (int j = 0; j < l_n; ++j) {
      sum += m[(static_cast<size_t>(i) * l_n) + j] * p[j];
    }
    ap[i] = sum;
  }
}

double SmyshlaevASleCgTaskALL::UpdateSolutionAndResidual(std::vector<double> &x, std::vector<double> &r,
                                                         const std::vector<double> &p, const std::vector<double> &ap,
                                                         double alpha, int start, int end, bool is_mpi) {
  double local_rs_new = 0.0;
#pragma omp parallel for default(none) shared(start, end, x, p, r, ap, alpha) reduction(+ : local_rs_new)
  for (int i = start; i < end; ++i) {
    x[i] += alpha * p[i];
    r[i] -= alpha * ap[i];
    local_rs_new += r[i] * r[i];
  }
  double global_rs_new = local_rs_new;
  if (is_mpi) {
    MPI_Allreduce(&local_rs_new, &global_rs_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  return global_rs_new;
}

void SmyshlaevASleCgTaskALL::SyncVectorP(std::vector<double> &p, int size, bool is_mpi) const {
  if (!is_mpi || n_ <= 0) {
    return;
  }
  std::vector<int> counts(size);
  std::vector<int> displs(size);
  for (int i = 0; i < size; ++i) {
    counts[i] = (n_ / size) + (i < (n_ % size) ? 1 : 0);
    displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
  }
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.data(), counts.data(), displs.data(), MPI_DOUBLE,
                 MPI_COMM_WORLD);
}

void SmyshlaevASleCgTaskALL::FinalGather(std::vector<double> &x, int start, int count, int size, bool is_mpi) {
  if (is_mpi && n_ > 0) {
    std::vector<int> counts(size);
    std::vector<int> displs(size);
    for (int i = 0; i < size; ++i) {
      counts[i] = (n_ / size) + (i < (n_ % size) ? 1 : 0);
      displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
    }
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      MPI_Gatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, x.data(), counts.data(), displs.data(), MPI_DOUBLE, 0,
                  MPI_COMM_WORLD);
      res_ = x;
    } else {
      MPI_Gatherv(x.data() + start, count, MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  } else {
    res_ = x;
  }
}

bool SmyshlaevASleCgTaskALL::RunImpl() {
  int rank = 0;
  int size = 1;
  bool is_mpi = ppc::util::IsUnderMpirun();
  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  std::vector<double> b_vec;
  if (rank == 0) {
    b_vec = GetInput().b;
  }
  DistributeInitialData(rank, is_mpi, b_vec);
  if (n_ == 0) {
    return true;
  }

  int l_n = n_;
  int start = (rank * (l_n / size)) + std::min(rank, l_n % size);
  int count = (l_n / size) + (rank < (l_n % size) ? 1 : 0);
  int end = start + count;

  std::vector<double> r = b_vec;
  std::vector<double> p = r;
  std::vector<double> x(l_n, 0.0);
  std::vector<double> ap(l_n, 0.0);
  omp_set_num_threads(ppc::util::GetNumThreads());

  double rs_old = ComputeDotProductAll(r, r, start, end, is_mpi);
  const double eps_sq = 1e-18;

  for (int iter = 0; iter < l_n * 2; ++iter) {
    if (rs_old < eps_sq) {
      break;
    }

    ComputeApAll(p, ap, start, end);
    double p_ap = ComputeDotProductAll(p, ap, start, end, is_mpi);
    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    double rs_new = UpdateSolutionAndResidual(x, r, p, ap, rs_old / p_ap, start, end, is_mpi);
    if (rs_new < eps_sq) {
      break;
    }

    double beta = rs_new / rs_old;
#pragma omp parallel for default(none) shared(start, end, p, r, beta)
    for (int i = start; i < end; ++i) {
      p[i] = r[i] + (beta * p[i]);
    }
    SyncVectorP(p, size, is_mpi);
    rs_old = rs_new;
  }
  FinalGather(x, start, count, size, is_mpi);
  return true;
}

bool SmyshlaevASleCgTaskALL::PostProcessingImpl() {
  int rank = 0;
  bool is_mpi = ppc::util::IsUnderMpirun();
  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Bcast(&n_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (n_ > 0) {
      if (rank != 0) {
        res_.resize(n_);
      }
      MPI_Bcast(res_.data(), n_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  }

  if (rank == 0 || is_mpi) {
    GetOutput() = res_;
  }
  return true;
}

}  // namespace smyshlaev_a_sle_cg_seq
