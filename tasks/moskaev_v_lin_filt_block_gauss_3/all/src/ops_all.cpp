// ops_all.cpp
#include "moskaev_v_lin_filt_block_gauss_3/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

namespace {

void CopyBlockWithHalo(const std::vector<uint8_t> &src, std::vector<uint8_t> &dst, int src_width, int src_height,
                       int channels, int block_x, int block_y, int block_w, int block_h, int padded_w) {
  for (int row = -1; row <= block_h; ++row) {
    for (int col = -1; col <= block_w; ++col) {
      int src_row = std::clamp(block_y + row, 0, src_height - 1);
      int src_col = std::clamp(block_x + col, 0, src_width - 1);
      int dst_row = row + 1;
      int dst_col = col + 1;
      for (int ch = 0; ch < channels; ++ch) {
        size_t src_idx = ((static_cast<size_t>(src_row) * src_width + src_col) * channels) + ch;
        size_t dst_idx = ((static_cast<size_t>(dst_row) * padded_w + dst_col) * channels) + ch;
        dst[dst_idx] = src[src_idx];
      }
    }
  }
}

void FilterBlock(const std::vector<uint8_t> &input_block, std::vector<uint8_t> &output_block, int block_w, int block_h,
                 int channels) {
  for (int row = 0; row < block_h; ++row) {
    for (int col = 0; col < block_w; ++col) {
      for (int ch = 0; ch < channels; ++ch) {
        float sum = 0.0F;
        for (int ky = -1; ky <= 1; ++ky) {
          for (int kx = -1; kx <= 1; ++kx) {
            int ny = row + 1 + ky;
            int nx = col + 1 + kx;
            size_t idx = ((static_cast<size_t>(ny) * (block_w + 2) + nx) * channels) + ch;
            int kidx = ((ky + 1) * 3) + (kx + 1);
            sum += static_cast<float>(input_block[idx]) * kGaussianKernel[kidx];
          }
        }
        size_t out_idx = ((static_cast<size_t>(row) * block_w + col) * channels) + ch;
        output_block[out_idx] = static_cast<uint8_t>(std::round(sum));
      }
    }
  }
}

void CopyBlockToOutput(const std::vector<uint8_t> &src_block, std::vector<uint8_t> &dst, int dst_width, int channels,
                       int block_x, int block_y, int block_w, int block_h) {
  for (int row = 0; row < block_h; ++row) {
    for (int col = 0; col < block_w; ++col) {
      for (int ch = 0; ch < channels; ++ch) {
        size_t src_idx = ((static_cast<size_t>(row) * block_w + col) * channels) + ch;
        size_t dst_idx = ((static_cast<size_t>(block_y + row) * dst_width + (block_x + col)) * channels) + ch;
        dst[dst_idx] = src_block[src_idx];
      }
    }
  }
}

}  // namespace

MoskaevVLinFiltBlockGauss3ALL::MoskaevVLinFiltBlockGauss3ALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
}

bool MoskaevVLinFiltBlockGauss3ALL::ValidationImpl() {
  const auto &input = GetInput();
  if (rank_ != 0) {
    return true;
  }
  const auto &data = std::get<4>(input);
  return !data.empty();
}

bool MoskaevVLinFiltBlockGauss3ALL::PreProcessingImpl() {
  return true;
}

bool MoskaevVLinFiltBlockGauss3ALL::PostProcessingImpl() {
  return !GetOutput().empty();
}

bool MoskaevVLinFiltBlockGauss3ALL::RunImpl() {
  int width = 0;
  int height = 0;
  int channels = 0;
  std::vector<uint8_t> image_data;

  if (rank_ == 0) {
    const auto &input = GetInput();
    width = std::get<0>(input);
    height = std::get<1>(input);
    channels = std::get<2>(input);
    image_data = std::get<4>(input);
  }

  MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&channels, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (width == 0 || height == 0) {
    return false;
  }

  size_t total_pixels = static_cast<size_t>(width) * height * channels;
  image_data.resize(total_pixels);
  MPI_Bcast(image_data.data(), static_cast<int>(total_pixels), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  std::vector<uint8_t> output(total_pixels, 0);

  int blocks_x = (width + block_size_ - 1) / block_size_;
  int blocks_y = (height + block_size_ - 1) / block_size_;
  int total_blocks = blocks_x * blocks_y;

  int blocks_per_proc = total_blocks / num_procs_;
  int remainder = total_blocks % num_procs_;
  int local_blocks = blocks_per_proc + (rank_ < remainder ? 1 : 0);

  std::vector<int> block_indices(total_blocks);
  for (int i = 0; i < total_blocks; ++i) {
    block_indices[i] = i;
  }

  std::vector<int> send_counts(num_procs_);
  std::vector<int> send_displs(num_procs_);
  int offset = 0;
  for (int proc = 0; proc < num_procs_; ++proc) {
    int cnt = blocks_per_proc + (proc < remainder ? 1 : 0);
    send_counts[proc] = cnt;
    send_displs[proc] = offset;
    offset += cnt;
  }

  std::vector<int> local_indices(local_blocks);
  MPI_Scatterv(block_indices.data(), send_counts.data(), send_displs.data(), MPI_INT, local_indices.data(),
               local_blocks, MPI_INT, 0, MPI_COMM_WORLD);

  for (int idx : local_indices) {
    int bx = idx % blocks_x;
    int by = idx / blocks_x;

    int block_x = bx * block_size_;
    int block_y = by * block_size_;
    int block_w = std::min(block_size_, width - block_x);
    int block_h = std::min(block_size_, height - block_y);

    int padded_w = block_w + 2;

    size_t input_block_size =
        static_cast<size_t>(padded_w) * static_cast<size_t>(block_h + 2) * static_cast<size_t>(channels);
    std::vector<uint8_t> input_block(input_block_size, 0);

    size_t output_block_size =
        static_cast<size_t>(block_w) * static_cast<size_t>(block_h) * static_cast<size_t>(channels);
    std::vector<uint8_t> output_block(output_block_size, 0);

    CopyBlockWithHalo(image_data, input_block, width, height, channels, block_x, block_y, block_w, block_h, padded_w);
    FilterBlock(input_block, output_block, block_w, block_h, channels);
    CopyBlockToOutput(output_block, output, width, channels, block_x, block_y, block_w, block_h);
  }

  std::vector<uint8_t> final_output;
  if (rank_ == 0) {
    final_output.resize(total_pixels);
  }

  MPI_Reduce(output.data(), final_output.data(), static_cast<int>(total_pixels), MPI_UNSIGNED_CHAR, MPI_BOR, 0,
             MPI_COMM_WORLD);

  if (rank_ == 0) {
    GetOutput() = std::move(final_output);
  }

  int output_size = static_cast<int>(GetOutput().size());
  MPI_Bcast(&output_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank_ != 0) {
    GetOutput().resize(output_size);
  }

  MPI_Bcast(GetOutput().data(), output_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  return true;
}

}  // namespace moskaev_v_lin_filt_block_gauss_3
