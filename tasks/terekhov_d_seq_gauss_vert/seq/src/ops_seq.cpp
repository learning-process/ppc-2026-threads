#include "terekhov_d_seq_gauss_vert/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "terekhov_d_seq_gauss_vert/common/include/common.hpp"

namespace terekhov_d_seq_gauss_vert {

TerekhovDGaussVertSEQ::TerekhovDGaussVertSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TerekhovDGaussVertSEQ::ValidationImpl() {
  const auto &input = GetInput();

  if (input.width <= 0 || input.height <= 0) {
    return false;
  }

  if (static_cast<int>(input.data.size()) != input.width * input.height) {
    return false;
  }

  return true;
}

bool TerekhovDGaussVertSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  width_ = input.width;
  height_ = input.height;

  GetOutput().width = width_;
  GetOutput().height = height_;
  GetOutput().data.resize(width_ * height_);

  int padded_width = width_ + 2;
  int padded_height = height_ + 2;
  padded_image_.resize(padded_width * padded_height);

  for (int y = 0; y < padded_height; ++y) {
    for (int x = 0; x < padded_width; ++x) {
      int src_x = x - 1;
      int src_y = y - 1;

      if (src_x < 0) {
        src_x = -src_x - 1;
      }
      if (src_x >= width_) {
        src_x = 2 * width_ - src_x - 1;
      }
      if (src_y < 0) {
        src_y = -src_y - 1;
      }
      if (src_y >= height_) {
        src_y = 2 * height_ - src_y - 1;
      }

      padded_image_[y * padded_width + x] = input.data[src_y * width_ + src_x];
    }
  }

  return true;
}

bool TerekhovDGaussVertSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (input.data.empty() || width_ <= 0 || height_ <= 0) {
    return false;
  }

  int padded_width = width_ + 2;

  const int num_bands = 4;
  int band_width = width_ / num_bands;
  if (band_width < 1) {
    band_width = 1;
  }

  for (int band = 0; band < num_bands; ++band) {
    int start_x = band * band_width;
    int end_x = (band == num_bands - 1) ? width_ : (band + 1) * band_width;

    for (int y = 0; y < height_; ++y) {
      for (int x = start_x; x < end_x; ++x) {
        int idx = y * width_ + x;

        float sum = 0.0f;

        for (int ky = -1; ky <= 1; ++ky) {
          for (int kx = -1; kx <= 1; ++kx) {
            int px = x + kx + 1;
            int py = y + ky + 1;

            int kernel_idx = (ky + 1) * 3 + (kx + 1);
            int pixel_value = padded_image_[py * padded_width + px];

            sum += pixel_value * kGaussKernel[kernel_idx];
          }
        }

        output.data[idx] = static_cast<int>(sum + 0.5f);
      }
    }
  }

  return true;
}

bool TerekhovDGaussVertSEQ::PostProcessingImpl() {
  return GetOutput().data.size() == static_cast<size_t>(GetOutput().width * GetOutput().height);
}

}  // namespace terekhov_d_seq_gauss_vert
