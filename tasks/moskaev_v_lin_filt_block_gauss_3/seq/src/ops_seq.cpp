#include "moskaev_v_lin_filt_block_gauss_3/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace moskaev_v_lin_filt_block_gauss_3 {

MoskaevVLinFiltBlockGauss3SEQ::MoskaevVLinFiltBlockGauss3SEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MoskaevVLinFiltBlockGauss3SEQ::ValidationImpl() {
  const auto &input = GetInput();
  return !std::get<4>(input).empty();
}

bool MoskaevVLinFiltBlockGauss3SEQ::PreProcessingImpl() {
  return true;
}

void MoskaevVLinFiltBlockGauss3SEQ::ApplyGaussianFilterToBlock(const std::vector<uint8_t> &input_block,
                                                               std::vector<uint8_t> &output_block, int block_width,
                                                               int block_height, int channels) {
  int inner_width = block_width - 2;
  int inner_height = block_height - 2;

  for (int y = 0; y < inner_height; ++y) {
    for (int x = 0; x < inner_width; ++x) {
      for (int c = 0; c < channels; ++c) {
        float sum = 0.0f;

        for (int ky = -1; ky <= 1; ++ky) {
          for (int kx = -1; kx <= 1; ++kx) {
            int ny = y + 1 + ky;
            int nx = x + 1 + kx;

            int idx = (ny * block_width + nx) * channels + c;
            sum += input_block[idx] * kGaussianKernel[(ky + 1) * 3 + (kx + 1)];
          }
        }

        int out_idx = (y * inner_width + x) * channels + c;
        output_block[out_idx] = static_cast<uint8_t>(std::round(sum));
      }
    }
  }
}

bool MoskaevVLinFiltBlockGauss3SEQ::RunImpl() {
  const auto &input = GetInput();

  int width = std::get<0>(input);
  int height = std::get<1>(input);
  int channels = std::get<2>(input);
  const auto &image_data = std::get<4>(input);

  if (image_data.empty()) {
    return false;
  }

  block_size_ = 64;

  GetOutput().resize(width * height * channels);

  for (int block_y = 0; block_y < height; block_y += block_size_) {
    for (int block_x = 0; block_x < width; block_x += block_size_) {
      int current_block_width = std::min(block_size_, width - block_x);
      int current_block_height = std::min(block_size_, height - block_y);

      int block_with_padding_width = current_block_width + 2;
      int block_with_padding_height = current_block_height + 2;

      std::vector<uint8_t> input_block(block_with_padding_width * block_with_padding_height * channels, 0);
      std::vector<uint8_t> output_block(current_block_width * current_block_height * channels, 0);

      for (int y = -1; y <= current_block_height; ++y) {
        for (int x = -1; x <= current_block_width; ++x) {
          int src_y = std::clamp(block_y + y, 0, height - 1);
          int src_x = std::clamp(block_x + x, 0, width - 1);
          int dst_y = y + 1;
          int dst_x = x + 1;

          for (int c = 0; c < channels; ++c) {
            int src_idx = (src_y * width + src_x) * channels + c;
            int dst_idx = (dst_y * block_with_padding_width + dst_x) * channels + c;
            input_block[dst_idx] = image_data[src_idx];
          }
        }
      }

      ApplyGaussianFilterToBlock(input_block, output_block, block_with_padding_width, block_with_padding_height,
                                 channels);

      for (int y = 0; y < current_block_height; ++y) {
        for (int x = 0; x < current_block_width; ++x) {
          for (int c = 0; c < channels; ++c) {
            int src_idx = (y * current_block_width + x) * channels + c;
            int dst_idx = ((block_y + y) * width + (block_x + x)) * channels + c;
            GetOutput()[dst_idx] = output_block[src_idx];
          }
        }
      }
    }
  }

  return true;
}

bool MoskaevVLinFiltBlockGauss3SEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace moskaev_v_lin_filt_block_gauss_3
