#include "moskaev_v_lin_filt_block_gauss_3/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

MoskaevVLinFiltBlockGauss3TBB::MoskaevVLinFiltBlockGauss3TBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MoskaevVLinFiltBlockGauss3TBB::ValidationImpl() {
  const auto &input = GetInput();
  return !std::get<4>(input).empty();
}

bool MoskaevVLinFiltBlockGauss3TBB::PreProcessingImpl() {
  return true;
}

namespace {

inline void ComputeFilteredPixel(const std::vector<uint8_t> &input_block, std::vector<uint8_t> &output_block,
                                 int block_width, int inner_width, int channels, int row, int col, int channel) {
  float sum = 0.0F;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      int ny = row + 1 + ky;
      int nx = col + 1 + kx;
      int idx = (((ny * block_width) + nx) * channels) + channel;
      sum += static_cast<float>(input_block[idx]) * kGaussianKernel[((ky + 1) * 3) + (kx + 1)];
    }
  }
  int out_idx = (((row * inner_width) + col) * channels) + channel;
  output_block[out_idx] = static_cast<uint8_t>(std::round(sum));
}

}  // namespace

void MoskaevVLinFiltBlockGauss3TBB::ApplyGaussianFilterToBlock(const std::vector<uint8_t> &input_block,
                                                               std::vector<uint8_t> &output_block, int block_width,
                                                               int block_height, int channels) {
  int inner_width = block_width - 2;
  int inner_height = block_height - 2;

  tbb::parallel_for(0, inner_height, [&](int row) {
    for (int col = 0; col < inner_width; ++col) {
      for (int channel = 0; channel < channels; ++channel) {
        ComputeFilteredPixel(input_block, output_block, block_width, inner_width, channels, row, col, channel);
      }
    }
  });
}

namespace {

void CopyBlockWithPadding(const std::vector<uint8_t> &source_image, std::vector<uint8_t> &padded_block, int width,
                          int height, int channels, int block_x, int block_y, int current_block_width,
                          int current_block_height, int block_with_padding_width) {
  int padded_height = current_block_height + 2;
  int padded_width = current_block_width + 2;

  tbb::parallel_for(0, padded_height, [&](int row) {
    for (int col = 0; col < padded_width; ++col) {
      int src_y = std::clamp(block_y + row - 1, 0, height - 1);
      int src_x = std::clamp(block_x + col - 1, 0, width - 1);
      for (int channel = 0; channel < channels; ++channel) {
        int src_idx = (((src_y * width) + src_x) * channels) + channel;
        int dst_idx = (((row * block_with_padding_width) + col) * channels) + channel;
        padded_block[dst_idx] = source_image[src_idx];
      }
    }
  });
}

void CopyProcessedBlockToOutput(const std::vector<uint8_t> &processed_block, std::vector<uint8_t> &output_image,
                                int width, int channels, int block_x, int block_y, int current_block_width,
                                int current_block_height) {
  tbb::parallel_for(0, current_block_height, [&](int row) {
    for (int col = 0; col < current_block_width; ++col) {
      for (int channel = 0; channel < channels; ++channel) {
        int src_idx = (((row * current_block_width) + col) * channels) + channel;
        int dst_idx = ((((block_y + row) * width) + (block_x + col)) * channels) + channel;
        output_image[dst_idx] = processed_block[src_idx];
      }
    }
  });
}

}  // namespace

bool MoskaevVLinFiltBlockGauss3TBB::RunImpl() {
  const auto &input = GetInput();

  int width = std::get<0>(input);
  int height = std::get<1>(input);
  int channels = std::get<2>(input);
  const auto &image_data = std::get<4>(input);

  if (image_data.empty()) {
    return false;
  }

  block_size_ = 64;
  int block_size = block_size_;

  GetOutput().resize(static_cast<size_t>(width) * static_cast<size_t>(height) * static_cast<size_t>(channels));

  int blocks_x = (width + block_size - 1) / block_size;
  int blocks_y = (height + block_size - 1) / block_size;

  tbb::parallel_for(0, blocks_y, [&](int by) {
    for (int bx = 0; bx < blocks_x; ++bx) {
      int block_x = bx * block_size;
      int block_y = by * block_size;
      int current_block_width = std::min(block_size, width - block_x);
      int current_block_height = std::min(block_size, height - block_y);
      int block_with_padding_width = current_block_width + 2;
      int block_with_padding_height = current_block_height + 2;

      std::vector<uint8_t> input_block(static_cast<size_t>(block_with_padding_width) *
                                           static_cast<size_t>(block_with_padding_height) *
                                           static_cast<size_t>(channels),
                                       0);
      std::vector<uint8_t> output_block(static_cast<size_t>(current_block_width) *
                                            static_cast<size_t>(current_block_height) * static_cast<size_t>(channels),
                                        0);

      CopyBlockWithPadding(image_data, input_block, width, height, channels, block_x, block_y, current_block_width,
                           current_block_height, block_with_padding_width);
      ApplyGaussianFilterToBlock(input_block, output_block, block_with_padding_width, block_with_padding_height,
                                 channels);
      CopyProcessedBlockToOutput(output_block, GetOutput(), width, channels, block_x, block_y, current_block_width,
                                 current_block_height);
    }
  });

  return true;
}

bool MoskaevVLinFiltBlockGauss3TBB::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace moskaev_v_lin_filt_block_gauss_3
