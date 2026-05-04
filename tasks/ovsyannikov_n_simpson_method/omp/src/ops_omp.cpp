#include "ovsyannikov_n_simpson_method/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"

namespace ovsyannikov_n_simpson_method {

double OvsyannikovNSimpsonMethodOMP::Function(double x, double y) {
  return x + y;
}

double OvsyannikovNSimpsonMethodOMP::GetCoeff(int i, int n) {
  if (i == 0 || i == n) {
    return 1.0;
  }
  return (i % 2 == 1) ? 4.0 : 2.0;
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
  const int nx_l = params_.nx;
  const int ny_l = params_.ny;
  const double ax_l = params_.ax;
  const double ay_l = params_.ay;
  const double hx_l = (params_.bx - params_.ax) / nx_l;
  const double hy_l = (params_.by - params_.ay) / ny_l;

  double total_sum = 0.0;

#pragma omp parallel for default(none) shared(nx_l, ny_l, ax_l, ay_l, hx_l, hy_l) reduction(+ : total_sum)
  for (int i = 0; i <= nx_l; ++i) {
    const double x = ax_l + (i * hx_l);
    const double coeff_x = GetCoeff(i, nx_l);
    double row_sum = 0.0;
    for (int j = 0; j <= ny_l; ++j) {
      const double y = ay_l + (j * hy_l);
      const double coeff_y = GetCoeff(j, ny_l);
      row_sum += coeff_y * Function(x, y);
    }
    total_sum += coeff_x * row_sum;
  }

  res_ = (hx_l * hy_l / 9.0) * total_sum;
  return true;
}

bool OvsyannikovNSimpsonMethodOMP::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace ovsyannikov_n_simpson_method
