#include "ovsyannikov_n_simpson_method/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <cmath>
#include <functional>

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"

namespace ovsyannikov_n_simpson_method {

double OvsyannikovNSimpsonMethodTBB::Function(double x, double y) {
  return x + y;
}

double OvsyannikovNSimpsonMethodTBB::GetCoeff(int i, int n) {
  if (i == 0 || i == n) {
    return 1.0;
  }
  return (i % 2 == 1) ? 4.0 : 2.0;
}

OvsyannikovNSimpsonMethodTBB::OvsyannikovNSimpsonMethodTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool OvsyannikovNSimpsonMethodTBB::ValidationImpl() {
  return GetInput().nx > 0 && GetInput().nx % 2 == 0 && GetInput().ny > 0 && GetInput().ny % 2 == 0;
}

bool OvsyannikovNSimpsonMethodTBB::PreProcessingImpl() {
  params_ = GetInput();
  res_ = 0.0;
  return true;
}

bool OvsyannikovNSimpsonMethodTBB::RunImpl() {
  const int nx_l = params_.nx;
  const int ny_l = params_.ny;
  const double ax_l = params_.ax;
  const double ay_l = params_.ay;
  const double hx = (params_.bx - params_.ax) / nx_l;
  const double hy = (params_.by - params_.ay) / ny_l;

  double total_sum = tbb::parallel_reduce(tbb::blocked_range<int>(0, nx_l + 1), 0.0,
                                          [&](const tbb::blocked_range<int> &r, double local_sum) {
    for (int i = r.begin(); i < r.end(); ++i) {
      const double x = ax_l + (i * hx);
      const double coeff_x = GetCoeff(i, nx_l);
      double row_sum = 0.0;
      for (int j = 0; j <= ny_l; ++j) {
        const double y = ay_l + (j * hy);
        const double coeff_y = GetCoeff(j, ny_l);
        row_sum += coeff_y * Function(x, y);
      }
      local_sum += coeff_x * row_sum;
    }
    return local_sum;
  }, std::plus<>());

  res_ = (hx * hy / 9.0) * total_sum;
  return true;
}

bool OvsyannikovNSimpsonMethodTBB::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}
}  // namespace ovsyannikov_n_simpson_method
