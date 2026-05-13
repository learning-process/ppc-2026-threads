#include "boltenkov_s_gaussian_kernel/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace boltenkov_s_gaussian_kernel {

BoltenkovSGaussianKernelSTL::BoltenkovSGaussianKernelSTL(const InType &in)
    : kernel_{{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}} {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::vector<int>>();
}

bool BoltenkovSGaussianKernelSTL::ValidationImpl() {
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

bool BoltenkovSGaussianKernelSTL::PreProcessingImpl() {
  GetOutput().resize(std::get<0>(GetInput()));
  for (std::size_t i = 0; i < std::get<0>(GetInput()); i++) {
    GetOutput()[i].resize(std::get<1>(GetInput()));
  }
  return true;
}

void BoltenkovSGaussianKernelSTL::CopyRowsToPaddedBuffer(const std::vector<std::vector<int>> &data,
                                                         std::vector<std::vector<int>> &tmp_data, std::size_t n,
                                                         unsigned int num_threads) {
  std::vector<std::thread> threads;
  std::size_t rows_per_thread = n / num_threads;
  std::size_t remainder = n % num_threads;
  std::size_t start = 1;

  for (unsigned int tid = 0; tid < num_threads; ++tid) {
    std::size_t count = rows_per_thread + (tid < remainder ? 1 : 0);
    std::size_t end = start + count;
    threads.emplace_back([&data, &tmp_data, start, end] {
      for (std::size_t i = start; i < end; ++i) {
        std::copy(data[i - 1].begin(), data[i - 1].end(), tmp_data[i].begin() + 1);
      }
    });
    start = end;
  }
  for (auto &th : threads) {
    th.join();
  }
}

void BoltenkovSGaussianKernelSTL::ApplyGaussianKernel(const std::vector<std::vector<int>> &tmp_data,
                                                      std::vector<std::vector<int>> &res, std::size_t n, std::size_t m,
                                                      unsigned int num_threads) {
  auto kernel = kernel_;
  int shift = shift_;
  std::vector<std::thread> threads;
  std::size_t rows_per_thread = n / num_threads;
  std::size_t remainder = n % num_threads;
  std::size_t start = 1;
  for (unsigned int tid = 0; tid < num_threads; ++tid) {
    std::size_t count = rows_per_thread + (tid < remainder ? 1 : 0);
    std::size_t end = start + count;
    threads.emplace_back([&tmp_data, &res, kernel, shift, start, end, m] {
      for (std::size_t i = start; i < end; ++i) {
        for (std::size_t j = 1; j <= m; ++j) {
          int sum = (tmp_data[i - 1][j - 1] * kernel[0][0]) + (tmp_data[i - 1][j] * kernel[0][1]) +
                    (tmp_data[i - 1][j + 1] * kernel[0][2]) + (tmp_data[i][j - 1] * kernel[1][0]) +
                    (tmp_data[i][j] * kernel[1][1]) + (tmp_data[i][j + 1] * kernel[1][2]) +
                    (tmp_data[i + 1][j - 1] * kernel[2][0]) + (tmp_data[i + 1][j] * kernel[2][1]) +
                    (tmp_data[i + 1][j + 1] * kernel[2][2]);
          res[i - 1][j - 1] = sum >> shift;
        }
      }
    });
    start = end;
  }
  for (auto &th : threads) {
    th.join();
  }
}

bool BoltenkovSGaussianKernelSTL::RunImpl() {
  std::size_t n = std::get<0>(GetInput());
  std::size_t m = std::get<1>(GetInput());
  std::vector<std::vector<int>> data = std::get<2>(GetInput());
  std::vector<std::vector<int>> tmp_data(n + 2, std::vector<int>(m + 2, 0));
  std::vector<std::vector<int>> &res = GetOutput();

  auto num_threads = static_cast<unsigned int>(ppc::util::GetNumThreads());
  if (num_threads == 0) {
    num_threads = 1;
  }

  CopyRowsToPaddedBuffer(data, tmp_data, n, num_threads);
  ApplyGaussianKernel(tmp_data, res, n, m, num_threads);

  return true;
}

bool BoltenkovSGaussianKernelSTL::PostProcessingImpl() {
  return true;
}

}  // namespace boltenkov_s_gaussian_kernel
