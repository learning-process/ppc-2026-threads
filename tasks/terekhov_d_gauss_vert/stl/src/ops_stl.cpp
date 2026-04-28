#include "terekhov_d_gauss_vert/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "terekhov_d_gauss_vert/common/include/common.hpp"
#include "util/include/util.hpp"

namespace terekhov_d_gauss_vert {

namespace {

inline void ProcessPixel(OutType &output, const std::vector<int> &padded_image, int padded_width, int width, int row,
                         int col) {
  size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col);
  float sum = 0.0F;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      int px = col + kx + 1;
      int py = row + ky + 1;
      int kernel_idx = ((ky + 1) * 3) + (kx + 1);
      size_t padded_idx = (static_cast<size_t>(py) * static_cast<size_t>(padded_width)) + static_cast<size_t>(px);
      sum += static_cast<float>(padded_image[padded_idx]) * kGaussKernel[static_cast<size_t>(kernel_idx)];
    }
  }
  output.data[idx] = static_cast<int>(std::lround(sum));
}

void ProcessBand(OutType &output, const std::vector<int> &padded_image, int padded_width, int width, int height,
                 int band, int band_width, int num_bands) {
  int start_x = band * band_width;
  int end_x = (band == num_bands - 1) ? width : ((band + 1) * band_width);
  for (int row = 0; row < height; ++row) {
    for (int col = start_x; col < end_x; ++col) {
      ProcessPixel(output, padded_image, padded_width, width, row, col);
    }
  }
}

OutType SolveSTL(const std::vector<int> &padded_image, int width, int height) {
  OutType output;
  output.width = width;
  output.height = height;
  output.data.resize(static_cast<size_t>(width) * static_cast<size_t>(height));

  const int padded_width = width + 2;
  const int num_threads = ppc::util::GetNumThreads();
  const int num_bands = std::min(num_threads, width);  // Исправлено!
  const int band_width = std::max(width / num_bands, 1);

  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_bands));

  for (int band = 0; band < num_bands; ++band) {
    threads.emplace_back([&output, &padded_image, padded_width, width, height, band, band_width, num_bands]() {
      ProcessBand(output, padded_image, padded_width, width, height, band, band_width, num_bands);
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  return output;
}

}  // namespace

TerekhovDGaussVertSTL::TerekhovDGaussVertSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TerekhovDGaussVertSTL::ValidationImpl() {
  const auto &input = GetInput();
  if (input.width <= 0 || input.height <= 0) {
    return false;
  }
  if (static_cast<int>(input.data.size()) != input.width * input.height) {
    return false;
  }
  return true;
}

bool TerekhovDGaussVertSTL::PreProcessingImpl() {
  const auto &input = GetInput();
  width_ = input.width;
  height_ = input.height;

  int padded_width = width_ + 2;
  int padded_height = height_ + 2;
  padded_image_.resize(static_cast<size_t>(padded_width) * static_cast<size_t>(padded_height));

  for (int row = 0; row < padded_height; ++row) {
    for (int col = 0; col < padded_width; ++col) {
      int src_x = col - 1;
      int src_y = row - 1;

      if (src_x < 0) {
        src_x = -src_x - 1;
      }
      if (src_x >= width_) {
        src_x = (2 * width_) - src_x - 1;
      }
      if (src_y < 0) {
        src_y = -src_y - 1;
      }
      if (src_y >= height_) {
        src_y = (2 * height_) - src_y - 1;
      }

      size_t padded_idx = (static_cast<size_t>(row) * static_cast<size_t>(padded_width)) + static_cast<size_t>(col);
      size_t src_idx = (static_cast<size_t>(src_y) * static_cast<size_t>(width_)) + static_cast<size_t>(src_x);
      padded_image_[padded_idx] = input.data[src_idx];
    }
  }
  return true;
}

bool TerekhovDGaussVertSTL::RunImpl() {
  const auto &input = GetInput();
  if (input.data.empty() || width_ <= 0 || height_ <= 0) {
    return false;
  }
  GetOutput() = SolveSTL(padded_image_, width_, height_);
  return true;
}

bool TerekhovDGaussVertSTL::PostProcessingImpl() {
  return GetOutput().data.size() == (static_cast<size_t>(GetOutput().width) * static_cast<size_t>(GetOutput().height));
}

}  // namespace terekhov_d_gauss_vert
