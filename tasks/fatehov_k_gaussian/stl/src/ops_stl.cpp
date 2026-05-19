#include "fatehov_k_gaussian/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "fatehov_k_gaussian/common/include/common.hpp"
#include "util/include/util.hpp"

namespace fatehov_k_gaussian {

namespace {

float ConvolvePixel(const Image &image, const std::vector<float> &kernel, int kernel_size, int half, int y_coord,
                    int x_coord, int c_coord) {
  const int w = static_cast<int>(image.width);
  const int h = static_cast<int>(image.height);
  const int ch = static_cast<int>(image.channels);
  float res = 0.0F;

  for (int ky = -half; ky <= half; ++ky) {
    for (int kx = -half; kx <= half; ++kx) {
      const int ny = std::clamp(y_coord + ky, 0, h - 1);
      const int nx = std::clamp(x_coord + kx, 0, w - 1);
      const float weight = kernel[((ky + half) * kernel_size) + (kx + half)];
      res += static_cast<float>(image.data[((ny * w + nx) * ch) + c_coord]) * weight;
    }
  }

  return res;
}

void ProcessPixel(const Image &src, Image &dst, const std::vector<float> &kernel, int kernel_size, int half,
                  int y_coord, int x_coord) {
  const int ch = static_cast<int>(src.channels);
  const int w = static_cast<int>(src.width);

  for (int c_coord = 0; c_coord < ch; ++c_coord) {
    const float res = ConvolvePixel(src, kernel, kernel_size, half, y_coord, x_coord, c_coord);
    dst.data[((y_coord * w + x_coord) * ch) + c_coord] = static_cast<uint8_t>(std::clamp(res, 0.0F, 255.0F));
  }
}

void ProcessRows(const Image &src, Image &dst, const std::vector<float> &kernel, int kernel_size, int half,
                 int row_begin, int row_end) {
  const int w = static_cast<int>(src.width);

  for (int y_coord = row_begin; y_coord < row_end; ++y_coord) {
    for (int x_coord = 0; x_coord < w; ++x_coord) {
      ProcessPixel(src, dst, kernel, kernel_size, half, y_coord, x_coord);
    }
  }
}

}  // namespace

FatehovKGaussianSTL::FatehovKGaussianSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FatehovKGaussianSTL::ValidationImpl() {
  const auto &input = GetInput();
  return input.image.width > 0 && input.image.height > 0 && input.image.channels > 0 && !input.image.data.empty() &&
         input.sigma > 0.0F;
}

bool FatehovKGaussianSTL::PreProcessingImpl() {
  const auto &input = GetInput();
  const float sigma = input.sigma;

  kernel_size_ = (2 * static_cast<int>(std::ceil(3.0F * sigma))) + 1;
  kernel_.resize(static_cast<std::size_t>(kernel_size_) * kernel_size_);

  const int half = kernel_size_ / 2;
  const float two_sigma_sq = 2.0F * sigma * sigma;
  float sum = 0.0F;

  for (int i = -half; i <= half; ++i) {
    for (int j = -half; j <= half; ++j) {
      const float val = std::exp(-(static_cast<float>((i * i) + (j * j))) / two_sigma_sq);
      kernel_[((i + half) * kernel_size_) + (j + half)] = val;
      sum += val;
    }
  }

  for (float &val : kernel_) {
    val /= sum;
  }

  GetOutput() = Image(input.image.width, input.image.height, input.image.channels);

  return true;
}

bool FatehovKGaussianSTL::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  const int h = static_cast<int>(input.image.height);
  const int half = kernel_size_ / 2;
  const int kernel_size = kernel_size_;
  const auto &kernel = kernel_;
  const int num_threads = ppc::util::GetNumThreads();

  std::vector<std::thread> threads(num_threads);
  const int rows_per_thread = h / num_threads;
  const int remainder = h % num_threads;

  int row_begin = 0;
  for (int t = 0; t < num_threads; ++t) {
    const int row_end = row_begin + rows_per_thread + (t < remainder ? 1 : 0);
    threads[t] = std::thread(ProcessRows, std::cref(input.image), std::ref(output), std::cref(kernel), kernel_size,
                             half, row_begin, row_end);
    row_begin = row_end;
  }

  for (auto &thread : threads) {
    thread.join();
  }

  return true;
}

bool FatehovKGaussianSTL::PostProcessingImpl() {
  return true;
}

}  // namespace fatehov_k_gaussian
