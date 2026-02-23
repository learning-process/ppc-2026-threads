#include "boltenkov_s_gaussian_kernel/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace boltenkov_s_gaussian_kernel {

BoltenkovSGaussianKernelSEQ::BoltenkovSGaussianKernelSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::vector<int>>();
  kernel = { {1, 2, 1}, {2, 4, 2}, {1, 2, 1} };
  shift = 4;
}

bool BoltenkovSGaussianKernelSEQ::ValidationImpl() {
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

bool BoltenkovSGaussianKernelSEQ::PreProcessingImpl() {
  GetOutput().resize(std::get<0>(GetInput()), std::vector<int>(std::get<1>(GetInput())));
  return true;
}

bool BoltenkovSGaussianKernelSEQ::RunImpl() {

  std::size_t n = std::get<0>(GetInput());
  std::size_t m = std::get<1>(GetInput());

  std::vector<std::vector<int>> data = std::get<2>(GetInput());
  std::vector<std::vector<int>> tmpData(n + 2, std::vector<int>(m + 2, 0));
  std::vector<std::vector<int>>& res = GetOutput();

  for (std::size_t i = 1; i <= n; i++) {
    std::copy(data[i - 1].begin(), data[i - 1].end(), tmpData[i].begin() + 1);
  }

  for (std::size_t i = 1; i <= n; i++) {
    for (std::size_t j = 1; j <= m; j++) {
      res[i - 1][j - 1] = tmpData[i - 1][j - 1] * kernel[0][0] + tmpData[i - 1][j] * kernel[0][1] + 
                          tmpData[i - 1][j + 1] * kernel[0][2] + tmpData[i][j - 1] * kernel[1][0] +
                          tmpData[i][j] * kernel[1][1] + tmpData[i][j + 1] * kernel[1][2] +
                          tmpData[i + 1][j - 1] * kernel[2][0] + tmpData[i + 1][j] * kernel[2][1] + 
                          tmpData[i + 1][j + 1] * kernel[2][2];
      res[i - 1][j - 1] >>= shift;
    }
  }

  return true;
}

bool BoltenkovSGaussianKernelSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace boltenkov_s_gaussian_kernel
