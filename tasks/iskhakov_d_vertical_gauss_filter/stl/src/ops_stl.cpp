#include "iskhakov_d_vertical_gauss_filter/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "iskhakov_d_vertical_gauss_filter/common/include/common.hpp"
#include "util/include/util.hpp"

namespace iskhakov_d_vertical_gauss_filter {

namespace {
const int kDivConst = 16;
const std::array<std::array<int, 3>, 3> kGaussKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t IskhakovDGetPixelMirrorSTL(const std::vector<uint8_t> &src, int col, int row, int width, int height) {
  if (col < 0) {
    col = -col - 1;
  } else if (col >= width) {
    col = (2 * width) - col - 1;
  }
  if (row < 0) {
    row = -row - 1;
  } else if (row >= height) {
    row = (2 * height) - row - 1;
  }
  return src[(row * width) + col];
}

}  // namespace

IskhakovDVerticalGaussFilterSTL::IskhakovDVerticalGaussFilterSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool IskhakovDVerticalGaussFilterSTL::ValidationImpl() {
  const auto &in = GetInput();
  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  if (in.data.size() != static_cast<size_t>(in.width) * static_cast<size_t>(in.height)) {
    return false;
  }
  return true;
}

bool IskhakovDVerticalGaussFilterSTL::PreProcessingImpl() {
  return true;
}

bool IskhakovDVerticalGaussFilterSTL::RunImpl() {
  const auto &in = GetInput();

  int width = in.width;
  int height = in.height;
  const std::vector<uint8_t> &matrix = in.data;
  std::vector<uint8_t> result(static_cast<size_t>(width) * static_cast<size_t>(height));

  const int num_threads = ppc::util::GetNumThreads();
  const int actual_threads = std::min(num_threads, width);
  std::vector<std::thread> threads;
  threads.reserve(actual_threads);

  const int cols_per_thread = width / actual_threads;
  const int remainder = width % actual_threads;
  int start_col = 0;

  for (int thread_id = 0; thread_id < actual_threads; ++thread_id) {
    int end_col = start_col + cols_per_thread + (thread_id < remainder ? 1 : 0);
    threads.emplace_back([&, start_col, end_col]() {
      for (int horizontal_band = start_col; horizontal_band < end_col; ++horizontal_band) {
        for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
          int sum = 0;

          sum += kGaussKernel[0][0] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band - 1, vertical_band - 1, width, height);
          sum += kGaussKernel[0][1] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band, vertical_band - 1, width, height);
          sum += kGaussKernel[0][2] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band + 1, vertical_band - 1, width, height);

          sum += kGaussKernel[1][0] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band - 1, vertical_band, width, height);
          sum += kGaussKernel[1][1] * IskhakovDGetPixelMirrorSTL(matrix, horizontal_band, vertical_band, width, height);
          sum += kGaussKernel[1][2] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band + 1, vertical_band, width, height);

          sum += kGaussKernel[2][0] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band - 1, vertical_band + 1, width, height);
          sum += kGaussKernel[2][1] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band, vertical_band + 1, width, height);
          sum += kGaussKernel[2][2] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band + 1, vertical_band + 1, width, height);

          result[(vertical_band * width) + horizontal_band] = static_cast<uint8_t>(sum / kDivConst);
        }
      }
    });
    start_col = end_col;
  }

  for (auto &th : threads) {
    th.join();
  }

  GetOutput().width = width;
  GetOutput().height = height;
  GetOutput().data = std::move(result);
  return true;
}

bool IskhakovDVerticalGaussFilterSTL::PostProcessingImpl() {
  return true;
}

}  // namespace iskhakov_d_vertical_gauss_filter
