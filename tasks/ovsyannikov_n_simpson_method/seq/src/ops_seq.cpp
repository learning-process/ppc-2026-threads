#include "ovsyannikov_n_simpson_method/seq/include/ops_seq.hpp"

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"

namespace ovsyannikov_n_simpson_method {

// Тестовая функция: f(x, y) = x + y
double OvsyannikovNSimpsonMethodSEQ::function(double x, double y) {
  return x + y;
}

OvsyannikovNSimpsonMethodSEQ::OvsyannikovNSimpsonMethodSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool OvsyannikovNSimpsonMethodSEQ::ValidationImpl() {
  return GetInput().nx > 0 && GetInput().nx % 2 == 0 && GetInput().ny > 0 && GetInput().ny % 2 == 0;
}

bool OvsyannikovNSimpsonMethodSEQ::PreProcessingImpl() {
  params_ = GetInput();
  res_ = 0;
  return true;
}

bool OvsyannikovNSimpsonMethodSEQ::RunImpl() {
  double hx = (params_.bx - params_.ax) / params_.nx;
  double hy = (params_.by - params_.ay) / params_.ny;
  double total_sum = 0.0;

  for (int i = 0; i <= params_.nx; ++i) {
    double x = params_.ax + (i * hx);
    double coeff_x = 2.0;
    if (i == 0 || i == params_.nx) {
      coeff_x = 1.0;
    } else if (i % 2 == 1) {
      coeff_x = 4.0;
    }

    for (int j = 0; j <= params_.ny; ++j) {
      double y = params_.ay + (j * hy);
      double coeff_y = 2.0;
      if (j == 0 || j == params_.ny) {
        coeff_y = 1.0;
      } else if (j % 2 == 1) {
        coeff_y = 4.0;
      }

      total_sum += coeff_x * coeff_y * function(x, y);
    }
  }

  res_ = (hx * hy / 9.0) * total_sum;
  return true;
}

bool OvsyannikovNSimpsonMethodSEQ::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace ovsyannikov_n_simpson_method
