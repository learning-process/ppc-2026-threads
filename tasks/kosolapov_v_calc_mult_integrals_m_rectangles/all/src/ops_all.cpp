#include "kosolapov_v_calc_mult_integrals_m_rectangles/all/include/ops_all.hpp"

#include <mpi.h>

#include <cmath>
#include <future>
#include <thread>
#include <tuple>
#include <vector>

#include "kosolapov_v_calc_mult_integrals_m_rectangles/common/include/common.hpp"

namespace kosolapov_v_calc_mult_integrals_m_rectangles {

KosolapovVCalcMultIntegralsMRectanglesALL::KosolapovVCalcMultIntegralsMRectanglesALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = InType(in);
  GetOutput() = 0.0;
}

bool KosolapovVCalcMultIntegralsMRectanglesALL::ValidationImpl() {
  int steps = std::get<0>(GetInput());
  int func_id = std::get<1>(GetInput());
  return steps > 0 && func_id >= 1 && func_id <= 4;
}

bool KosolapovVCalcMultIntegralsMRectanglesALL::PreProcessingImpl() {
  return true;
}

bool KosolapovVCalcMultIntegralsMRectanglesALL::RunImpl() {
  int steps = std::get<0>(GetInput());
  int func_id = std::get<1>(GetInput());
  std::tuple<double, double, double, double> temp = GetBounds(func_id);
  double a = std::get<0>(temp);
  double b = std::get<1>(temp);
  double c = std::get<2>(temp);
  double d = std::get<3>(temp);
  double integral = RectanglesIntegral(func_id, steps, a, b, c, d);
  GetOutput() = integral;
  return true;
}

bool KosolapovVCalcMultIntegralsMRectanglesALL::PostProcessingImpl() {
  return true;
}

double KosolapovVCalcMultIntegralsMRectanglesALL::Function1(double x, double y) {
  // f(x,y) = x^2 + y^2
  return (x * x) + (y * y);
}
double KosolapovVCalcMultIntegralsMRectanglesALL::Function2(double x, double y) {
  // f(x,y) = sin(x) * cos(y)
  return std::sin(x) * std::cos(y);
}
double KosolapovVCalcMultIntegralsMRectanglesALL::Function3(double x, double y) {
  // f(x,y) = exp(-(x^2 + y^2))
  return std::exp(-((x * x) + (y * y)));
}
double KosolapovVCalcMultIntegralsMRectanglesALL::Function4(double x, double y) {
  // f(x,y) = sin(x + y)
  return std::sin(x + y);
}
double KosolapovVCalcMultIntegralsMRectanglesALL::CallFunction(int func_id, double x, double y) {
  switch (func_id) {
    case 1:
      return Function1(x, y);
    case 2:
      return Function2(x, y);
    case 3:
      return Function3(x, y);
    case 4:
      return Function4(x, y);
    default:
      return Function1(x, y);
  }
}
std::tuple<double, double, double, double> KosolapovVCalcMultIntegralsMRectanglesALL::GetBounds(int func_id) {
  switch (func_id) {
    case 1:
      return {0.0, 1.0, 0.0, 1.0};
    case 2:
      return {0.0, kPi, 0.0, kPi / 2.0};
    case 3:
      return {-1.0, 1.0, -1.0, 1.0};
    case 4:
      return {0.0, kPi, 0.0, kPi};
    default:
      return {0.0, 1.0, 0.0, 1.0};
  }
}
double KosolapovVCalcMultIntegralsMRectanglesALL::RectanglesIntegral(int func_id, int steps, double a, double b,
                                                                     double c, double d) {
  double hx = (b - a) / steps;
  double hy = (d - c) / steps;
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int rows_per_proc = steps / size;
  int remainder = steps % size;
  int start_i = rank * rows_per_proc + std::min(rank, remainder);
  int end_i = start_i + rows_per_proc + (rank < remainder ? 1 : 0);
  int my_rows = end_i - start_i;
  double local_sum = 0.0;
  if (my_rows > 0) {
    unsigned int num_threads = std::max(1U, std::thread::hardware_concurrency());
    if (num_threads > static_cast<unsigned int>(my_rows)) {
      num_threads = static_cast<unsigned int>(my_rows);
    }
    std::vector<std::future<double>> futures;
    futures.reserve(num_threads);
    int rows_per_thread = my_rows / num_threads;
    int thread_remainder = my_rows % num_threads;
    int offset = 0;
    for (unsigned int t = 0; t < num_threads; ++t) {
      int t_start = offset;
      int t_end = t_start + rows_per_thread + (t < static_cast<unsigned int>(thread_remainder) ? 1 : 0);
      offset = t_end;
      futures.push_back(std::async(std::launch::async, [=]() -> double {
        double part_sum = 0.0;
        for (int i = start_i + t_start; i < start_i + t_end; ++i) {
          double x = a + (i + 0.5) * hx;
          for (int j = 0; j < steps; ++j) {
            double y = c + (j + 0.5) * hy;
            part_sum += CallFunction(func_id, x, y);
          }
        }
        return part_sum;
      }));
    }
    for (auto &fut : futures) {
      local_sum += fut.get();
    }
  }
  double global_sum = 0.0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&global_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return global_sum * hx * hy;
}

}  // namespace kosolapov_v_calc_mult_integrals_m_rectangles
