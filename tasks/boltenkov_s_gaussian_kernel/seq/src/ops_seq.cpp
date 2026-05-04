#include "boltenkov_s_gaussian_kernel/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"

namespace boltenkov_s_gaussian_kernel {

BoltenkovSGaussianKernelSEQ::BoltenkovSGaussianKernelSEQ(const InType &in)
    : kernel_{{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}} {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::vector<int>>();
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
  GetOutput().resize(std::get<0>(GetInput()));
  for (std::size_t i = 0; i < std::get<0>(GetInput()); i++) {
    GetOutput()[i].resize(std::get<1>(GetInput()));
  }
  return true;
}

bool BoltenkovSGaussianKernelSEQ::RunImpl() {
  std::size_t n = std::get<0>(GetInput());
  std::size_t m = std::get<1>(GetInput());

  std::vector<std::vector<int>> data = std::get<2>(GetInput());
  std::vector<std::vector<int>> tmp_data(n + 2, std::vector<int>(m + 2, 0));
  std::vector<std::vector<int>> &res = GetOutput();

  for (std::size_t i = 1; i <= n; i++) {
    std::copy(data[i - 1].begin(), data[i - 1].end(), tmp_data[i].begin() + 1);
  }

  for (std::size_t i = 1; i <= n; i++) {
    for (std::size_t j = 1; j <= m; j++) {
      res[i - 1][j - 1] = (tmp_data[i - 1][j - 1] * kernel_[0][0]) + (tmp_data[i - 1][j] * kernel_[0][1]) +
                          (tmp_data[i - 1][j + 1] * kernel_[0][2]) + (tmp_data[i][j - 1] * kernel_[1][0]) +
                          (tmp_data[i][j] * kernel_[1][1]) + (tmp_data[i][j + 1] * kernel_[1][2]) +
                          (tmp_data[i + 1][j - 1] * kernel_[2][0]) + (tmp_data[i + 1][j] * kernel_[2][1]) +
                          (tmp_data[i + 1][j + 1] * kernel_[2][2]);
      res[i - 1][j - 1] >>= shift_;
    }
  }

  return true;
}

bool BoltenkovSGaussianKernelSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace boltenkov_s_gaussian_kernel
