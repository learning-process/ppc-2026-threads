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
  const auto &input = GetInput();
  int width = std::get<0>(input);
  int height = std::get<1>(input);
  int channels = std::get<2>(input);
  const auto &data = std::get<4>(input);
  return !data.empty() && static_cast<size_t>(width) * height * channels == data.size();
}

bool MoskaevVLinFiltBlockGauss3ALL::PreProcessingImpl() {
  return true;
}

namespace {

constexpr int BLOCK_SIZE = 256;

inline uint8_t ComputeFilteredPixel(const std::vector<uint8_t> &image, int width, int height, int channels, int row,
                                    int col, int channel) {
  float sum = 0.0F;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      int ny = std::clamp(row + ky, 0, height - 1);
      int nx = std::clamp(col + kx, 0, width - 1);
      int idx = (((ny * width) + nx) * channels) + channel;
      sum += static_cast<float>(image[idx]) * kGaussianKernel[((ky + 1) * 3) + (kx + 1)];
    }
  }
  return static_cast<uint8_t>(std::round(sum));
}

void ProcessBlock(const std::vector<uint8_t> &input, std::vector<uint8_t> &output, int width, int height, int channels,
                  int start_row, int block_x, int block_y, int block_width, int block_height) {
  for (int row = 0; row < block_height; ++row) {
    int global_row = start_row + block_y + row;
    for (int col = 0; col < block_width; ++col) {
      int global_col = block_x + col;
      for (int ch = 0; ch < channels; ++ch) {
        int output_row = block_y + row;
        int output_col = block_x + col;
        int output_idx = ((output_row * width) + output_col) * channels + ch;
        output[output_idx] = ComputeFilteredPixel(input, width, height, channels, global_row, global_col, ch);
      }
    }
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

  MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&channels, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int total_pixels = width * height * channels;
  std::vector<uint8_t> local_image(total_pixels);
  if (rank_ == 0) {
    local_image = image_data;
  }
  MPI_Bcast(local_image.data(), total_pixels, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  int rows_per_proc = height / num_procs_;
  int remainder = height % num_procs_;

  std::vector<int> send_counts(num_procs_);
  std::vector<int> displs(num_procs_);
  int offset = 0;
  for (int proc = 0; proc < num_procs_; ++proc) {
    int rows = rows_per_proc + (proc < remainder ? 1 : 0);
    send_counts[proc] = rows * width * channels;
    displs[proc] = offset;
    offset += send_counts[proc];
  }

  int local_rows = rows_per_proc + (rank_ < remainder ? 1 : 0);
  int local_size = local_rows * width * channels;
  std::vector<uint8_t> local_output(local_size, 0);

  int start_row = 0;
  for (int proc = 0; proc < rank_; ++proc) {
    start_row += rows_per_proc + (proc < remainder ? 1 : 0);
  }

  int blocks_x = (width + BLOCK_SIZE - 1) / BLOCK_SIZE;
  int blocks_y = (local_rows + BLOCK_SIZE - 1) / BLOCK_SIZE;

#pragma omp parallel for collapse(2) schedule(dynamic) default(none) \
    shared(local_image, local_output, width, height, channels, local_rows, start_row, blocks_x, blocks_y)
  for (int by = 0; by < blocks_y; ++by) {
    for (int bx = 0; bx < blocks_x; ++bx) {
      int block_x = bx * BLOCK_SIZE;
      int block_y = by * BLOCK_SIZE;
      int block_width = std::min(BLOCK_SIZE, width - block_x);
      int block_height = std::min(BLOCK_SIZE, local_rows - block_y);

      ProcessBlock(local_image, local_output, width, height, channels, start_row, block_x, block_y, block_width,
                   block_height);
    }
  }

  if (rank_ == 0) {
    GetOutput().resize(total_pixels);
  }
  MPI_Gatherv(local_output.data(), local_size, MPI_UNSIGNED_CHAR, rank_ == 0 ? GetOutput().data() : nullptr,
              send_counts.data(), displs.data(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  return true;
}

bool MoskaevVLinFiltBlockGauss3ALL::PostProcessingImpl() {
  return rank_ == 0 ? !GetOutput().empty() : true;
}

}  // namespace moskaev_v_lin_filt_block_gauss_3
