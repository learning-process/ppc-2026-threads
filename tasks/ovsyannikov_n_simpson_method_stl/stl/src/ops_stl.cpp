#include "ovsyannikov_n_simpson_method_stl/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <functional>
#include <numeric>
#include <vector>

#include "ovsyannikov_n_simpson_method_stl/common/include/common.hpp"

namespace ovsyannikov_n_simpson_method_stl {

double OvsyannikovNSimpsonMethodSTL::Function(double x, double y) {
  return x + y;
}

double OvsyannikovNSimpsonMethodSTL::GetCoeff(int i, int n) {
  if (i == 0 || i == n) {
    return 1.0;
  }
  return (i % 2 == 1) ? 4.0 : 2.0;
}

OvsyannikovNSimpsonMethodSTL::OvsyannikovNSimpsonMethodSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool OvsyannikovNSimpsonMethodSTL::ValidationImpl() {
  return GetInput().nx > 0 && GetInput().nx % 2 == 0 && GetInput().ny > 0 && GetInput().ny % 2 == 0;
}

bool OvsyannikovNSimpsonMethodSTL::PreProcessingImpl() {
  params_ = GetInput();
  res_ = 0.0;
  return true;
}

bool OvsyannikovNSimpsonMethodSTL::RunImpl() {
  const int nx_l = params_.nx;
  const int ny_l = params_.ny;
  const double ax_l = params_.ax;
  const double ay_l = params_.ay;
  const double hx = (params_.bx - params_.ax) / nx_l;
  const double hy = (params_.by - params_.ay) / ny_l;

  std::vector<int> indices(nx_l + 1);
  std::iota(indices.begin(), indices.end(), 0);

  double total_sum =
      std::transform_reduce(std::execution::par, indices.begin(), indices.end(), 0.0, std::plus<>(), [&](int i) {
    const double x = ax_l + (i * hx);
    const double coeff_x = GetCoeff(i, nx_l);
    double row_sum = 0.0;
    for (int j = 0; j <= ny_l; ++j) {
      const double y = ay_l + (j * hy);
      const double coeff_y = GetCoeff(j, ny_l);
      row_sum += coeff_y * Function(x, y);
    }
    return coeff_x * row_sum;
  });

  res_ = (hx * hy / 9.0) * total_sum;
  return true;
}

bool OvsyannikovNSimpsonMethodSTL::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}
}  // namespace ovsyannikov_n_simpson_method_stl
