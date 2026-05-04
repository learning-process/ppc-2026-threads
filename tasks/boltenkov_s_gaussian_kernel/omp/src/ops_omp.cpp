#include "boltenkov_s_gaussian_kernel/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace boltenkov_s_gaussian_kernel {

BoltenkovSGaussianKernelOMP::BoltenkovSGaussianKernelOMP(const InType &in)
    : kernel_{{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}} {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::vector<int>>();
}

bool BoltenkovSGaussianKernelOMP::ValidationImpl() {
  std::size_t n = std::get<0>(GetInput());
  std::size_t m = std::get<1>(GetInput());
  if (std::get<2>(GetInput()).size() != n) {
    return false;
  }
  for (std::size_t i = 0; i < n; i++) {
    if (std::get<2>(GetInput())[i].size() != m) {
      return false;
    }
  }
  return true;
}

bool BoltenkovSGaussianKernelOMP::PreProcessingImpl() {
  GetOutput().resize(std::get<0>(GetInput()));
  for (std::size_t i = 0; i < std::get<0>(GetInput()); i++) {
    GetOutput()[i].resize(std::get<1>(GetInput()));
  }
  return true;
}

bool BoltenkovSGaussianKernelOMP::RunImpl() {
  std::size_t n = std::get<0>(GetInput());
  std::size_t m = std::get<1>(GetInput());

  std::vector<std::vector<int>> data = std::get<2>(GetInput());
  std::vector<std::vector<int>> tmp_data(n + 2, std::vector<int>(m + 2, 0));
  std::vector<std::vector<int>> &res = GetOutput();

#pragma omp parallel for num_threads(ppc::util::GetNumThreads()) default(none) shared(tmp_data, data) firstprivate(n)
  for (std::size_t i = 1; i <= n; i++) {
    std::copy(data[i - 1].begin(), data[i - 1].end(), tmp_data[i].begin() + 1);
  }

  auto kernel = kernel_;
  int shift = shift_;

#pragma omp parallel for num_threads(ppc::util::GetNumThreads()) default(none) shared(tmp_data, res) \
    firstprivate(n, m, kernel, shift)
  for (std::size_t i = 1; i <= n; i++) {
    for (std::size_t j = 1; j <= m; j++) {
      res[i - 1][j - 1] = (tmp_data[i - 1][j - 1] * kernel[0][0]) + (tmp_data[i - 1][j] * kernel[0][1]) +
                          (tmp_data[i - 1][j + 1] * kernel[0][2]) + (tmp_data[i][j - 1] * kernel[1][0]) +
                          (tmp_data[i][j] * kernel[1][1]) + (tmp_data[i][j + 1] * kernel[1][2]) +
                          (tmp_data[i + 1][j - 1] * kernel[2][0]) + (tmp_data[i + 1][j] * kernel[2][1]) +
                          (tmp_data[i + 1][j + 1] * kernel[2][2]);
      res[i - 1][j - 1] >>= shift;
    }
  }

  return true;
}

bool BoltenkovSGaussianKernelOMP::PostProcessingImpl() {
  return true;
}

}  // namespace boltenkov_s_gaussian_kernel
