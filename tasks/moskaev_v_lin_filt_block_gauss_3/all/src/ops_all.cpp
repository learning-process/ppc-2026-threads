#include "moskaev_v_lin_filt_block_gauss_3/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

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
                                int current_block_height) {
  for (int row = 0; row < current_block_height; ++row) {
    for (int col = 0; col < current_block_width; ++col) {
      for (int channel = 0; channel < channels; ++channel) {
        int src_idx = (((row * current_block_width) + col) * channels) + channel;
        int dst_idx = ((((block_y + row) * width) + (block_x + col)) * channels) + channel;
        output_image[dst_idx] = processed_block[src_idx];
      }
    }
  }
}

int ComputeStartBlock(int rank, int blocks_per_proc, int remainder) {
  return (rank < remainder) ? rank * (blocks_per_proc + 1) : (rank * blocks_per_proc) + remainder;
}

int ComputeEndBlock(int rank, int blocks_per_proc, int remainder) {
  return ComputeStartBlock(rank, blocks_per_proc, remainder) + blocks_per_proc + (rank < remainder ? 1 : 0);
}

void ProcessBlocks(const std::vector<uint8_t> &local_image, std::vector<uint8_t> &output, int width, int height,
                   int channels, int block_size, int start_block_y, int end_block_y) {
  for (int by = start_block_y; by < end_block_y; ++by) {
    for (int bx = 0; bx < (width + block_size - 1) / block_size; ++bx) {
      int block_x = bx * block_size;
      int block_y = by * block_size;
      int current_block_width = std::min(block_size, width - block_x);
      int current_block_height = std::min(block_size, height - block_y);
      int pw = current_block_width + 2;
      int ph = current_block_height + 2;

      std::vector<uint8_t> in_block(static_cast<size_t>(pw) * ph * channels, 0);
      std::vector<uint8_t> out_block(static_cast<size_t>(current_block_width) * current_block_height * channels, 0);

      CopyBlockWithPadding(local_image, in_block, width, height, channels, block_x, block_y, current_block_width,
                           current_block_height, pw);
      MoskaevVLinFiltBlockGauss3ALL::ApplyGaussianFilterToBlock(in_block, out_block, pw, ph, channels);
      CopyProcessedBlockToOutput(out_block, output, width, channels, block_x, block_y, current_block_width,
                                 current_block_height);
    }
  }
}

void GatherResults(std::vector<uint8_t> &output, int width, int channels, int block_size, int num_procs, int rank,
                   int blocks_per_proc, int remainder, int start_block_y, int end_block_y) {
  if (rank == 0) {
    std::vector<uint8_t> recv_buf;
    for (int proc = 1; proc < num_procs; ++proc) {
      int proc_start = ComputeStartBlock(proc, blocks_per_proc, remainder);
      int proc_end = ComputeEndBlock(proc, blocks_per_proc, remainder);
      int recv_count = width * (proc_end - proc_start) * block_size * channels;
      recv_buf.resize(recv_count);
      MPI_Recv(recv_buf.data(), recv_count, MPI_UNSIGNED_CHAR, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int offset = proc_start * block_size * width * channels;
      std::ranges::copy(recv_buf, output.begin() + offset);
    }
  } else {
    int send_count = width * (end_block_y - start_block_y) * block_size * channels;
    MPI_Send(output.data(), send_count, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
  }
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

  block_size_ = 64;
  int block_size = block_size_;

  if (rank_ == 0) {
    GetOutput().resize(static_cast<size_t>(width) * height * channels);
  }

  int total_size = static_cast<int>(image_data.size());
  MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<uint8_t> local_image(total_size);
  if (rank_ == 0) {
    local_image = image_data;
  }
  MPI_Bcast(local_image.data(), total_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  int blocks_y = (height + block_size - 1) / block_size;
  int blocks_per_proc = blocks_y / num_procs_;
  int remainder = blocks_y % num_procs_;

  int start_block_y = ComputeStartBlock(rank_, blocks_per_proc, remainder);
  int end_block_y = ComputeEndBlock(rank_, blocks_per_proc, remainder);

  if (rank_ != 0) {
    GetOutput().resize(static_cast<size_t>(width) * (end_block_y - start_block_y) * block_size * channels);
  }

  ProcessBlocks(local_image, GetOutput(), width, height, channels, block_size, start_block_y, end_block_y);
  GatherResults(GetOutput(), width, channels, block_size, num_procs_, rank_, blocks_per_proc, remainder, start_block_y,
                end_block_y);

  return true;
}

bool MoskaevVLinFiltBlockGauss3ALL::PostProcessingImpl() {
  return rank_ == 0 ? !GetOutput().empty() : true;
}

}  // namespace moskaev_v_lin_filt_block_gauss_3
