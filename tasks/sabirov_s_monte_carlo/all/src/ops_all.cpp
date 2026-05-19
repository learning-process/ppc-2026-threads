#include "sabirov_s_monte_carlo/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "sabirov_s_monte_carlo/common/include/common.hpp"
#include "util/include/util.hpp"

namespace sabirov_s_monte_carlo {

namespace {

double EvalLinear(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x;
  }
  return s;
}

double EvalSumCubes(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x * x * x;
  }
  return s;
}

double EvalCosProduct(const std::vector<double> &point) {
  double p = 1.0;
  for (double x : point) {
    p *= std::cos(x);
  }
  return p;
}

double EvalExpNeg(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x;
  }
  return std::exp(-s);
}

double EvalMixedPoly(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += (x * x) + x;
  }
  return s;
}

double EvalSinSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += std::sin(x);
  }
  return s;
}

double EvalSqrtSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += std::sqrt(x);
  }
  return s;
}

double EvalQuarticSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x * x * x * x;
  }
  return s;
}

double EvaluateAt(FuncType func_type, const std::vector<double> &point) {
  switch (func_type) {
    case FuncType::kLinear:
      return EvalLinear(point);
    case FuncType::kSumCubes:
      return EvalSumCubes(point);
    case FuncType::kCosProduct:
      return EvalCosProduct(point);
    case FuncType::kExpNeg:
      return EvalExpNeg(point);
    case FuncType::kMixedPoly:
      return EvalMixedPoly(point);
    case FuncType::kSinSum:
      return EvalSinSum(point);
    case FuncType::kSqrtSum:
      return EvalSqrtSum(point);
    case FuncType::kQuarticSum:
      return EvalQuarticSum(point);
    default:
      return 0.0;
  }
}

void FillScatterCounts(int num_samples, int size, std::vector<int> *counts) {
  const int base = num_samples / size;
  const int rem = num_samples % size;
  for (int rk = 0; rk < size; ++rk) {
    (*counts)[static_cast<size_t>(rk)] = base + (rk < rem ? 1 : 0);
  }
}

void SendParamsToNonRoot(int size, const std::vector<double> &lower, const std::vector<double> &upper, int num_samples,
                         FuncType func_type) {
  const int dims_send = static_cast<int>(lower.size());
  for (int dest = 1; dest < size; ++dest) {
    std::array<int, 3> head = {dims_send, num_samples, static_cast<int>(func_type)};
    MPI_Send(head.data(), 3, MPI_INT, dest, 0, MPI_COMM_WORLD);
    MPI_Send(lower.data(), dims_send, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
    MPI_Send(upper.data(), dims_send, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
  }
}

void RecvParamsFromRoot(std::vector<double> *lower, std::vector<double> *upper, int *num_samples, FuncType *func_type) {
  std::array<int, 3> head{};
  MPI_Recv(head.data(), 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  const int rdims = head[0];
  *num_samples = head[1];
  *func_type = static_cast<FuncType>(head[2]);
  lower->resize(static_cast<size_t>(rdims));
  upper->resize(static_cast<size_t>(rdims));
  MPI_Recv(lower->data(), rdims, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(upper->data(), rdims, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

double HyperrectangleVolume(const std::vector<double> &lower, const std::vector<double> &upper) {
  double volume = 1.0;
  const int dims = static_cast<int>(lower.size());
  for (int vi = 0; vi < dims; ++vi) {
    volume *= (upper[static_cast<size_t>(vi)] - lower[static_cast<size_t>(vi)]);
  }
  return volume;
}

double ParallelLocalSum(const std::vector<double> &lower, const std::vector<double> &upper, int local_n,
                        FuncType ftype) {
  const int dims = static_cast<int>(lower.size());
  const std::vector<double> *const plower = &lower;
  const std::vector<double> *const pupper = &upper;
  const FuncType copy_ftype = ftype;
  const int ln = local_n;
  double sum = 0.0;

#pragma omp parallel default(none) shared(plower, pupper, copy_ftype, dims, ln) reduction(+ : sum)
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<std::uniform_real_distribution<double>> dists;
    dists.reserve(static_cast<size_t>(dims));
    for (int wi = 0; wi < dims; ++wi) {
      dists.emplace_back((*plower)[static_cast<size_t>(wi)], (*pupper)[static_cast<size_t>(wi)]);
    }
    std::vector<double> point(static_cast<size_t>(dims));

#pragma omp for schedule(static)
    for (int si = 0; si < ln; ++si) {
      for (int sj = 0; sj < dims; ++sj) {
        point[static_cast<size_t>(sj)] = dists[static_cast<size_t>(sj)](gen);
      }
      sum += EvaluateAt(copy_ftype, point);
    }
  }
  return sum;
}

}  // namespace

SabirovSMonteCarloALL::SabirovSMonteCarloALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SabirovSMonteCarloALL::ValidationImpl() {
  const auto &in = GetInput();
  if (in.lower.size() != in.upper.size() || in.lower.empty()) {
    return false;
  }
  if (in.num_samples <= 0) {
    return false;
  }
  for (size_t i = 0; i < in.lower.size(); ++i) {
    if (in.lower[i] >= in.upper[i]) {
      return false;
    }
  }
  if (in.func_type < FuncType::kLinear || in.func_type > FuncType::kQuarticSum) {
    return false;
  }
  constexpr size_t kMaxDimensions = 10;
  return in.lower.size() <= kMaxDimensions;
}

bool SabirovSMonteCarloALL::PreProcessingImpl() {
  const auto &in = GetInput();
  lower_ = in.lower;
  upper_ = in.upper;
  num_samples_ = in.num_samples;
  func_type_ = in.func_type;
  GetOutput() = 0.0;
  return true;
}

bool SabirovSMonteCarloALL::RunImpl() {
  omp_set_num_threads(std::max(1, ppc::util::GetNumThreads()));

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> counts(static_cast<size_t>(size));
  if (rank == 0) {
    FillScatterCounts(num_samples_, size, &counts);
  }

  int local_n = 0;
  std::vector<int> scounts(static_cast<size_t>(size), 1);
  std::vector<int> displs(static_cast<size_t>(size));
  for (int di = 0; di < size; ++di) {
    displs[static_cast<size_t>(di)] = di;
  }

  const int *scatter_send = (rank == 0) ? counts.data() : nullptr;
  MPI_Scatterv(scatter_send, scounts.data(), displs.data(), MPI_INT, &local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int global_n = num_samples_;
  if (rank == 0) {
    SendParamsToNonRoot(size, lower_, upper_, num_samples_, func_type_);
  } else {
    RecvParamsFromRoot(&lower_, &upper_, &num_samples_, &func_type_);
    global_n = num_samples_;
  }

  const double volume = HyperrectangleVolume(lower_, upper_);
  double sum = ParallelLocalSum(lower_, upper_, local_n, func_type_);

  double global_sum = 0.0;
  MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  GetOutput() = volume * global_sum / static_cast<double>(global_n);
  return true;
}

bool SabirovSMonteCarloALL::PostProcessingImpl() {
  return true;
}

}  // namespace sabirov_s_monte_carlo
