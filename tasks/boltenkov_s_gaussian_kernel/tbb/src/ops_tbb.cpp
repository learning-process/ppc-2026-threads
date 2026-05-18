#include "boltenkov_s_gaussian_kernel/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace boltenkov_s_gaussian_kernel {

BoltenkovSGaussianKernelTBB::BoltenkovSGaussianKernelTBB(const InType &in)
    : kernel_{{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}} {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::vector<int>>();
}

bool BoltenkovSGaussianKernelTBB::ValidationImpl() {
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

bool BoltenkovSGaussianKernelTBB::PreProcessingImpl() {
  GetOutput().resize(std::get<0>(GetInput()));
  for (std::size_t i = 0; i < std::get<0>(GetInput()); i++) {
    GetOutput()[i].resize(std::get<1>(GetInput()));
  }
  return true;
}

bool BoltenkovSGaussianKernelTBB::RunImpl() {
  std::size_t n = std::get<0>(GetInput());
  std::size_t m = std::get<1>(GetInput());
  std::vector<std::vector<int>> data = std::get<2>(GetInput());
  std::vector<std::vector<int>> tmp_data(n + 2, std::vector<int>(m + 2, 0));
  std::vector<std::vector<int>> &res = GetOutput();

  tbb::task_arena arena(ppc::util::GetNumThreads());

  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<std::size_t>(1, n + 1), [&](const tbb::blocked_range<std::size_t> &r) {
      for (std::size_t i = r.begin(); i != r.end(); ++i) {
        std::copy(data[i - 1].begin(), data[i - 1].end(), tmp_data[i].begin() + 1);
      }
    });
  });

  auto kernel = kernel_;
  int shift = shift_;

  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<std::size_t>(1, n + 1), [&](const tbb::blocked_range<std::size_t> &r) {
      for (std::size_t i = r.begin(); i != r.end(); ++i) {
        for (std::size_t j = 1; j <= m; ++j) {
          res[i - 1][j - 1] = (tmp_data[i - 1][j - 1] * kernel[0][0]) + (tmp_data[i - 1][j] * kernel[0][1]) +
                              (tmp_data[i - 1][j + 1] * kernel[0][2]) + (tmp_data[i][j - 1] * kernel[1][0]) +
                              (tmp_data[i][j] * kernel[1][1]) + (tmp_data[i][j + 1] * kernel[1][2]) +
                              (tmp_data[i + 1][j - 1] * kernel[2][0]) + (tmp_data[i + 1][j] * kernel[2][1]) +
                              (tmp_data[i + 1][j + 1] * kernel[2][2]);
          res[i - 1][j - 1] >>= shift;
        }
      }
    });
  });

  return true;
}

bool BoltenkovSGaussianKernelTBB::PostProcessingImpl() {
  return true;
}

}  // namespace boltenkov_s_gaussian_kernel
