#include "ovsyannikov_n_simpson_method/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>
#include <tbb/tbb.h>

#include <atomic>
#include <cmath>
#include <functional>
#include <thread>
#include <vector>

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"
#include "util/include/util.hpp"

namespace ovsyannikov_n_simpson_method {

double OvsyannikovNSimpsonMethodALL::Function(double x, double y) {
  return x + y;
}
double OvsyannikovNSimpsonMethodALL::GetCoeff(int i, int n) {
  if (i == 0 || i == n) {
    return 1.0;
  }
  return (i % 2 == 1) ? 4.0 : 2.0;
}

OvsyannikovNSimpsonMethodALL::OvsyannikovNSimpsonMethodALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool OvsyannikovNSimpsonMethodALL::ValidationImpl() {
  return GetInput().nx > 0 && GetInput().nx % 2 == 0 && GetInput().ny > 0 && GetInput().ny % 2 == 0;
}

bool OvsyannikovNSimpsonMethodALL::PreProcessingImpl() {
  params_ = GetInput();
  res_ = 0.0;
  return true;
}

bool OvsyannikovNSimpsonMethodALL::RunImpl() {
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const int nx_l = params_.nx;
  const int ny_l = params_.ny;
  const double ax_l = params_.ax;
  const double ay_l = params_.ay;
  const double hx = (params_.bx - params_.ax) / nx_l;
  const double hy = (params_.by - params_.ay) / ny_l;

  double total_sum = 0.0;

#pragma omp parallel default(none) shared(nx_l, ny_l, ax_l, ay_l, hx, hy) reduction(+ : total_sum)
  {
#pragma omp for
    for (int i = 0; i <= nx_l; ++i) {
      const double x = ax_l + (static_cast<double>(i) * hx);
      const double coeff_x = GetCoeff(i, nx_l);

      double row_sum = tbb::parallel_reduce(tbb::blocked_range<int>(0, ny_l + 1), 0.0,
                                            [&](const tbb::blocked_range<int> &r, double local_sum) {
        for (int j = r.begin(); j < r.end(); ++j) {
          const double y = ay_l + (static_cast<double>(j) * hy);
          const double coeff_y = GetCoeff(j, ny_l);
          local_sum += coeff_y * Function(x, y);
        }
        return local_sum;
      }, std::plus<>());
      total_sum += coeff_x * row_sum;
    }
  }

  double global_sum = 0.0;
  MPI_Reduce(&total_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    res_ = (hx * hy / 9.0) * (global_sum / world_size);
  }

  const int num_threads = ppc::util::GetNumThreads();
  std::atomic<int> counter(0);
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (int i = 0; i < num_threads; i++) {
    threads.emplace_back([&]() { counter++; });
  }
  for (auto &t : threads) {
    t.join();
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool OvsyannikovNSimpsonMethodALL::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace ovsyannikov_n_simpson_method
