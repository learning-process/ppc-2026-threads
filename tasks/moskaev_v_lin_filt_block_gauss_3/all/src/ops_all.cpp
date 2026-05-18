#include "moskaev_v_lin_filt_block_gauss_3/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

MoskaevVLinFiltBlockGauss3ALL::MoskaevVLinFiltBlockGauss3ALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
}

bool MoskaevVLinFiltBlockGauss3ALL::ValidationImpl() {
  return !std::get<4>(GetInput()).empty();
}

bool MoskaevVLinFiltBlockGauss3ALL::PreProcessingImpl() {
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

void MoskaevVLinFiltBlockGauss3ALL::ApplyGaussianFilterToBlock(const std::vector<uint8_t> &input_block,
                                                               std::vector<uint8_t> &output_block, int block_width,
                                                               int block_height, int channels) {
  int inner_width = block_width - 2;
  int inner_height = block_height - 2;

#pragma omp parallel for collapse(3) schedule(static) default(none) \
    shared(input_block, output_block, inner_width, inner_height, channels, block_width)
  for (int row = 0; row < inner_height; ++row) {
    for (int col = 0; col < inner_width; ++col) {
      for (int channel = 0; channel < channels; ++channel) {
        ComputeFilteredPixel(input_block, output_block, block_width, inner_width, channels, row, col, channel);
      }
    }
  }
}

namespace {

void CopyBlockWithPadding(const std::vector<uint8_t> &source_image, std::vector<uint8_t> &padded_block, int width,
                          int height, int channels, int block_x, int block_y, int current_block_width,
                          int current_block_height, int block_with_padding_width) {
  for (int row = -1; row <= current_block_height; ++row) {
    for (int col = -1; col <= current_block_width; ++col) {
      int src_y = std::clamp(block_y + row, 0, height - 1);
      int src_x = std::clamp(block_x + col, 0, width - 1);
      int dst_y = row + 1;
      int dst_x = col + 1;
      for (int channel = 0; channel < channels; ++channel) {
        int src_idx = (((src_y * width) + src_x) * channels) + channel;
        int dst_idx = (((dst_y * block_with_padding_width) + dst_x) * channels) + channel;
        padded_block[dst_idx] = source_image[src_idx];
      }
    }
  }
}

void CopyProcessedBlockToOutput(const std::vector<uint8_t> &processed_block, std::vector<uint8_t> &output_image,
                                int width, int channels, int block_x, int block_y, int current_block_width,
                                int current_block_height, int output_start_row) {
  for (int row = 0; row < current_block_height; ++row) {
    for (int col = 0; col < current_block_width; ++col) {
      for (int channel = 0; channel < channels; ++channel) {
        int src_idx = (((row * current_block_width) + col) * channels) + channel;
        int dst_row = output_start_row + block_y + row;
        int dst_idx = ((dst_row * width) + (block_x + col)) * channels + channel;
        output_image[dst_idx] = processed_block[src_idx];
      }
    }
  }
}

int ComputeStartBlock(int rank, int blocks_per_proc, int remainder) {
  if (rank < remainder) {
    return rank * (blocks_per_proc + 1);
  }
  return rank * blocks_per_proc + remainder;
}

int ComputeEndBlock(int rank, int blocks_per_proc, int remainder) {
  int start = ComputeStartBlock(rank, blocks_per_proc, remainder);
  int extra = (rank < remainder) ? 1 : 0;
  return start + blocks_per_proc + extra;
}

void ProcessSingleBlock(const std::vector<uint8_t> &image, std::vector<uint8_t> &output, int width, int height,
                        int channels, int block_size, int bx, int by, int output_start_row) {
  int block_x = bx * block_size;
  int block_y = by * block_size;
  int current_block_width = std::min(block_size, width - block_x);
  int current_block_height = std::min(block_size, height - block_y);
  int pw = current_block_width + 2;
  int ph = current_block_height + 2;

  std::vector<uint8_t> in_block(static_cast<size_t>(pw) * ph * channels, 0);
  std::vector<uint8_t> out_block(static_cast<size_t>(current_block_width) * current_block_height * channels, 0);

  CopyBlockWithPadding(image, in_block, width, height, channels, block_x, block_y, current_block_width,
                       current_block_height, pw);
  MoskaevVLinFiltBlockGauss3ALL::ApplyGaussianFilterToBlock(in_block, out_block, pw, ph, channels);
  CopyProcessedBlockToOutput(out_block, output, width, channels, block_x, block_y, current_block_width,
                             current_block_height, output_start_row);
}

void ProcessBlockRange(const std::vector<uint8_t> &image, std::vector<uint8_t> &output, int width, int height,
                       int channels, int block_size, int start_block_y, int end_block_y, int output_start_row) {
  int blocks_x = (width + block_size - 1) / block_size;
  for (int by = start_block_y; by < end_block_y; ++by) {
    for (int bx = 0; bx < blocks_x; ++bx) {
      ProcessSingleBlock(image, output, width, height, channels, block_size, bx, by, output_start_row);
    }
  }
}

void BroadcastImage(int rank, const std::vector<uint8_t> &image_data, std::vector<uint8_t> &local_image) {
  int total_size = static_cast<int>(image_data.size());
  MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  local_image.resize(total_size);
  if (rank == 0) {
    local_image = image_data;
  }
  MPI_Bcast(local_image.data(), total_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
}

void SendResults(int dest, const std::vector<uint8_t> &data, size_t size) {
  MPI_Send(const_cast<uint8_t *>(data.data()), static_cast<int>(size), MPI_UNSIGNED_CHAR, dest, 0, MPI_COMM_WORLD);
}

void ReceiveResults(int src, std::vector<uint8_t> &buffer, size_t size) {
  buffer.resize(size);
  MPI_Recv(buffer.data(), static_cast<int>(size), MPI_UNSIGNED_CHAR, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

}  // namespace

bool MoskaevVLinFiltBlockGauss3ALL::RunImpl() {
  const auto &input = GetInput();
  int width = std::get<0>(input);
  int height = std::get<1>(input);
  int channels = std::get<2>(input);
  const auto &image_data = std::get<4>(input);

  if (image_data.empty()) {
    return false;
  }

  block_size_ = 256;
  int block_size = block_size_;

  std::vector<uint8_t> local_image;
  BroadcastImage(rank_, image_data, local_image);

  int blocks_y = (height + block_size - 1) / block_size;

  if (blocks_y <= rank_) {
    if (rank_ == 0) {
      GetOutput().resize(static_cast<size_t>(width) * height * channels);
      int blocks_x = (width + block_size - 1) / block_size;
      for (int by = 0; by < blocks_y; ++by) {
        for (int bx = 0; bx < blocks_x; ++bx) {
          ProcessSingleBlock(local_image, GetOutput(), width, height, channels, block_size, bx, by, 0);
        }
      }
    }
    return true;
  }

  int blocks_per_proc = blocks_y / num_procs_;
  int remainder = blocks_y % num_procs_;

  int start_block_y = ComputeStartBlock(rank_, blocks_per_proc, remainder);
  int end_block_y = ComputeEndBlock(rank_, blocks_per_proc, remainder);

  int start_row = start_block_y * block_size;
  int end_row = std::min(end_block_y * block_size, height);
  int local_height = end_row - start_row;

  size_t total_size = static_cast<size_t>(width) * height * channels;
  size_t local_size = static_cast<size_t>(width) * local_height * channels;

  if (rank_ == 0) {
    GetOutput().resize(total_size);
    ProcessBlockRange(local_image, GetOutput(), width, height, channels, block_size, start_block_y, end_block_y, 0);

    for (int proc = 1; proc < num_procs_; ++proc) {
      int p_start = ComputeStartBlock(proc, blocks_per_proc, remainder);
      int p_end = ComputeEndBlock(proc, blocks_per_proc, remainder);
      int p_start_row = p_start * block_size;
      int p_end_row = std::min(p_end * block_size, height);
      size_t p_size = static_cast<size_t>(width) * (p_end_row - p_start_row) * channels;

      if (p_size == 0) {
        continue;
      }

      std::vector<uint8_t> proc_data;
      ReceiveResults(proc, proc_data, p_size);

      int dst_offset = p_start_row * width * channels;
      std::copy(proc_data.begin(), proc_data.end(), GetOutput().begin() + dst_offset);
    }
  } else {
    if (local_size == 0) {
      return true;
    }

    std::vector<uint8_t> local_output(local_size);
    ProcessBlockRange(local_image, local_output, width, height, channels, block_size, start_block_y, end_block_y,
                      start_row);
    SendResults(0, local_output, local_size);
  }

  return true;
}

bool MoskaevVLinFiltBlockGauss3ALL::PostProcessingImpl() {
  if (rank_ == 0) {
    return !GetOutput().empty();
  }
  return true;
}

}  // namespace moskaev_v_lin_filt_block_gauss_3
