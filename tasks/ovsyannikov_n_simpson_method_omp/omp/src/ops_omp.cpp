#include "ovsyannikov_n_simpson_method_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>

namespace ovsyannikov_n_simpson_method_omp {
double OvsyannikovNSimpsonMethodOMP::Function(double x, double y) {
  return x + y;
}

OvsyannikovNSimpsonMethodOMP::OvsyannikovNSimpsonMethodOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool OvsyannikovNSimpsonMethodOMP::ValidationImpl() {
  return GetInput().nx > 0 && GetInput().nx % 2 == 0 && GetInput().ny > 0 && GetInput().ny % 2 == 0;
}

bool OvsyannikovNSimpsonMethodOMP::PreProcessingImpl() {
  params_ = GetInput();
  res_ = 0.0;
  return true;
}

bool OvsyannikovNSimpsonMethodOMP::RunImpl() {
  double hx = (params_.bx - params_.ax) / params_.nx;
  double hy = (params_.by - params_.ay) / params_.ny;
  double total_sum = 0.0;

#pragma omp parallel for reduction(+ : total_sum)
  for (int i = 0; i <= params_.nx; ++i) {
    double x = params_.ax + i * hx;
    double coeff_x = (i == 0 || i == params_.nx) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
    double row_sum = 0.0;
    for (int j = 0; j <= params_.ny; ++j) {
      double y = params_.ay + j * hy;
      double coeff_y = (j == 0 || j == params_.ny) ? 1.0 : (j % 2 == 1 ? 4.0 : 2.0);
      row_sum += coeff_y * Function(x, y);
    }
    total_sum += coeff_x * row_sum;
  }
  res_ = (hx * hy / 9.0) * total_sum;
  return true;
}

bool OvsyannikovNSimpsonMethodOMP::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}
}  // namespace ovsyannikov_n_simpson_method_omp
