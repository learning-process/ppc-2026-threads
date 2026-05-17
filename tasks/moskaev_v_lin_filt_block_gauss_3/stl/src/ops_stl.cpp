#include "moskaev_v_lin_filt_block_gauss_3/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <numeric>
#include <vector>

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

MoskaevVLinFiltBlockGauss3STL::MoskaevVLinFiltBlockGauss3STL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MoskaevVLinFiltBlockGauss3STL::ValidationImpl() {
  const auto &input = GetInput();
  return !std::get<4>(input).empty();
}

bool MoskaevVLinFiltBlockGauss3STL::PreProcessingImpl() {
  return true;
}

void MoskaevVLinFiltBlockGauss3STL::ApplyGaussianFilterToBlock(const std::vector<uint8_t> &input_block,
                                                               std::vector<uint8_t> &output_block, int block_width,
                                                               int block_height, int channels) {
  int inner_width = block_width - 2;
  int inner_height = block_height - 2;

  int total_pixels = inner_height * inner_width;

  std::vector<int> indices(total_pixels);
  std::ranges::iota(indices, 0);

  std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [&](int linear_idx) {
    int row = linear_idx / inner_width;
    int col = linear_idx % inner_width;

    for (int channel = 0; channel < channels; ++channel) {
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
  });
}

namespace {
void CopyBlockWithPadding(const std::vector<uint8_t> &source_image, std::vector<uint8_t> &padded_block, int width,
                          int height, int channels, int block_x, int block_y, int current_block_width,
                          int current_block_height, int block_with_padding_width) {
  int total_rows = current_block_height + 2;
  int total_cols = current_block_width + 2;

  std::vector<int> linear_indices(static_cast<size_t>(total_rows) * static_cast<size_t>(total_cols));
  std::ranges::iota(linear_indices, 0);

  std::for_each(std::execution::par_unseq, linear_indices.begin(), linear_indices.end(), [&](int linear_idx) {
    int dst_row = linear_idx / total_cols;
    int dst_col = linear_idx % total_cols;
    int src_row_unclamped = block_y + dst_row - 1;
    int src_col_unclamped = block_x + dst_col - 1;

    int src_y = std::clamp(src_row_unclamped, 0, height - 1);
    int src_x = std::clamp(src_col_unclamped, 0, width - 1);

    for (int channel = 0; channel < channels; ++channel) {
      int src_idx = (((src_y * width) + src_x) * channels) + channel;
      int dst_idx = (((dst_row * block_with_padding_width) + dst_col) * channels) + channel;
      padded_block[dst_idx] = source_image[src_idx];
    }
  });
}

void CopyProcessedBlockToOutput(const std::vector<uint8_t> &processed_block, std::vector<uint8_t> &output_image,
                                int width, int channels, int block_x, int block_y, int current_block_width,
                                int current_block_height) {
  std::vector<int> linear_indices(static_cast<size_t>(current_block_height) * static_cast<size_t>(current_block_width));
  std::ranges::iota(linear_indices, 0);

  std::for_each(std::execution::par_unseq, linear_indices.begin(), linear_indices.end(), [&](int linear_idx) {
    int row = linear_idx / current_block_width;
    int col = linear_idx % current_block_width;

    for (int channel = 0; channel < channels; ++channel) {
      int src_idx = (((row * current_block_width) + col) * channels) + channel;
      int dst_idx = ((((block_y + row) * width) + (block_x + col)) * channels) + channel;
      output_image[dst_idx] = processed_block[src_idx];
    }
  });
}
}  // namespace

bool MoskaevVLinFiltBlockGauss3STL::RunImpl() {
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

  std::vector<int> block_indices(static_cast<size_t>(blocks_y) * static_cast<size_t>(blocks_x));
  std::ranges::iota(block_indices, 0);

  std::for_each(std::execution::par_unseq, block_indices.begin(), block_indices.end(), [&](int linear_idx) {
    int by = linear_idx / blocks_x;
    int bx = linear_idx % blocks_x;

    int block_x = bx * block_size;
    int block_y = by * block_size;

    int current_block_width = std::min(block_size, width - block_x);
    int current_block_height = std::min(block_size, height - block_y);

    int block_with_padding_width = current_block_width + 2;
    int block_with_padding_height = current_block_height + 2;

    std::vector<uint8_t> input_block(static_cast<size_t>(block_with_padding_width) *
                                         static_cast<size_t>(block_with_padding_height) * static_cast<size_t>(channels),
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
  });

  return true;
}

bool MoskaevVLinFiltBlockGauss3STL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace moskaev_v_lin_filt_block_gauss_3
