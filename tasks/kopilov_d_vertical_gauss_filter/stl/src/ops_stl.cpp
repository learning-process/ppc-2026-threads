#include "kopilov_d_vertical_gauss_filter/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kopilov_d_vertical_gauss_filter {

namespace {
const int kDivisor = 16;
const std::array<std::array<int, 3>, 3> kGaussKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t GetPixelMirroredSTL(const std::vector<uint8_t> &src, int column, int row, int width, int height) {
  int new_col = column;
  int new_row = row;
  if (new_col < 0) {
    new_col = -new_col - 1;
  } else if (new_col >= width) {
    new_col = (2 * width) - new_col - 1;
  }
  if (new_row < 0) {
    new_row = -new_row - 1;
  } else if (new_row >= height) {
    new_row = (2 * height) - new_row - 1;
  }
  auto idx = (static_cast<size_t>(new_row) * static_cast<size_t>(width)) + static_cast<size_t>(new_col);
  return src[idx];
}

uint8_t ApplyGaussKernel(const std::vector<uint8_t> &src, int column, int row, int width, int height) {
  int pixel_sum = 0;
  for (size_t kernel_y = 0; kernel_y < 3; ++kernel_y) {
    for (size_t kernel_x = 0; kernel_x < 3; ++kernel_x) {
      int current_column = column + static_cast<int>(kernel_x) - 1;
      int current_row = row + static_cast<int>(kernel_y) - 1;
      pixel_sum +=
          kGaussKernel.at(kernel_y).at(kernel_x) * GetPixelMirroredSTL(src, current_column, current_row, width, height);
    }
  }
  return static_cast<uint8_t>(pixel_sum / kDivisor);
}
}  // namespace

KopilovDVerticalGaussFilterSTL::KopilovDVerticalGaussFilterSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool KopilovDVerticalGaussFilterSTL::ValidationImpl() {
  const auto &in = GetInput();
  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  return in.data.size() == static_cast<size_t>(in.width) * static_cast<size_t>(in.height);
}

bool KopilovDVerticalGaussFilterSTL::PreProcessingImpl() {
  return true;
}

bool KopilovDVerticalGaussFilterSTL::RunImpl() {
  const auto &in = GetInput();
  int width = in.width;
  int height = in.height;
  const std::vector<uint8_t> &matrix = in.data;
  std::vector<uint8_t> result(static_cast<size_t>(width) * static_cast<size_t>(height));

  const int num_threads = ppc::util::GetNumThreads();
  const int actual_threads_count = std::min(num_threads, width);

  std::vector<std::thread> threads;
  threads.reserve(actual_threads_count);

  const int cols_per_thread = width / actual_threads_count;
  const int remainder_cols = width % actual_threads_count;
  int start_column = 0;

  for (int thread_idx = 0; thread_idx < actual_threads_count; ++thread_idx) {
    int end_column = start_column + cols_per_thread + (thread_idx < remainder_cols ? 1 : 0);

    threads.emplace_back([&, start_column, end_column]() {
      for (int col_idx = start_column; col_idx < end_column; ++col_idx) {
        for (int row_idx = 0; row_idx < height; ++row_idx) {
          auto out_idx = (static_cast<size_t>(row_idx) * static_cast<size_t>(width)) + static_cast<size_t>(col_idx);
          result[out_idx] = ApplyGaussKernel(matrix, col_idx, row_idx, width, height);
        }
      }
    });
    start_column = end_column;
  }

  for (auto &thread_obj : threads) {
    thread_obj.join();
  }

  GetOutput().width = width;
  GetOutput().height = height;
  GetOutput().data = std::move(result);
  return true;
}

bool KopilovDVerticalGaussFilterSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kopilov_d_vertical_gauss_filter
