#include "romanov_a_gauss_block/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "romanov_a_gauss_block/common/include/common.hpp"
#include "util/include/util.hpp"

namespace romanov_a_gauss_block {

RomanovAGaussBlockALL::RomanovAGaussBlockALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  }
  GetOutput() = std::vector<uint8_t>();
}

bool RomanovAGaussBlockALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }
  return std::get<0>(GetInput()) * std::get<1>(GetInput()) * 3 ==
         static_cast<int>(std::get<2>(GetInput()).size());
}

bool RomanovAGaussBlockALL::PreProcessingImpl() {
  return true;
}

namespace {

constexpr int kBlockSize = 32;

int ApplyKernel(const std::vector<uint8_t> &img, int row, int col, int channel,
                int width, int buffer_height, int halo_top,
                const std::array<std::array<int, 3>, 3> &kernel) {
  int sum = 0;
  for (size_t kr = 0; kr < 3; ++kr) {
    for (size_t kc = 0; kc < 3; ++kc) {
      int nr_local = row + static_cast<int>(kr) - 1;
      int nc = col + static_cast<int>(kc) - 1;
      int buffer_row = nr_local + halo_top;
      if (buffer_row >= 0 && buffer_row < buffer_height && nc >= 0 && nc < width) {
        size_t idx = (((static_cast<size_t>(buffer_row) * width) + nc) * 3) + channel;
        sum += (static_cast<int>(img[idx]) * kernel.at(kr).at(kc));
      }
    }
  }
  return sum;
}

void ProcessFullBlock(const std::vector<uint8_t> &input, std::vector<uint8_t> &output,
                      int width, int buffer_height, int halo_top, int start_row, int start_col) {
  static constexpr std::array<std::array<int, 3>, 3> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

  for (int row = start_row; row < start_row + kBlockSize; ++row) {
    for (int col = start_col; col < start_col + kBlockSize; ++col) {
      for (int channel = 0; channel < 3; ++channel) {
        int sum = ApplyKernel(input, row, col, channel, width, buffer_height, halo_top, kKernel);
        int result_value = (sum + 8) / 16;
        result_value = std::clamp(result_value, 0, 255);
        auto idx = ((static_cast<size_t>(row) * width + col) * 3) + channel;
        output[idx] = static_cast<uint8_t>(result_value);
      }
    }
  }
}

void ProcessPartBlock(const std::vector<uint8_t> &input, std::vector<uint8_t> &output,
                      int width, int local_rows, int buffer_height, int halo_top,
                      int start_row, int start_col) {
  static constexpr std::array<std::array<int, 3>, 3> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

  const int end_row = std::min(local_rows, start_row + kBlockSize);
  const int end_col = std::min(width, start_col + kBlockSize);

  for (int row = start_row; row < end_row; ++row) {
    for (int col = start_col; col < end_col; ++col) {
      for (int channel = 0; channel < 3; ++channel) {
        int sum = ApplyKernel(input, row, col, channel, width, buffer_height, halo_top, kKernel);
        int result_value = (sum + 8) / 16;
        result_value = std::clamp(result_value, 0, 255);
        auto idx = ((static_cast<size_t>(row) * width + col) * 3) + channel;
        output[idx] = static_cast<uint8_t>(result_value);
      }
    }
  }
}

}  // namespace

bool RomanovAGaussBlockALL::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  std::array<int, 2> dims{};
  if (rank == 0) {
    dims[0] = std::get<0>(GetInput());
    dims[1] = std::get<1>(GetInput());
  }
  MPI_Bcast(dims.data(), 2, MPI_INT, 0, MPI_COMM_WORLD);
  const int width = dims[0];
  const int height = dims[1];

  const int total_block_rows = height / kBlockSize;
  const int height_remainder = height % kBlockSize;
  const int num_col_blocks = width / kBlockSize;
  const bool width_has_remainder = (width % kBlockSize) != 0;

  std::vector<int> block_rows_per_proc(world_size);
  std::vector<int> block_row_displs(world_size);
  {
    const int base_blocks = total_block_rows / world_size;
    const int extra_blocks = total_block_rows % world_size;
    int block_offset = 0;
    for (int p = 0; p < world_size; ++p) {
      block_rows_per_proc[p] = base_blocks + (p < extra_blocks ? 1 : 0);
      block_row_displs[p] = block_offset;
      block_offset += block_rows_per_proc[p];
    }
  }

  std::vector<int> rows_per_proc(world_size);
  std::vector<int> row_displs(world_size);
  {
    int pixel_offset = 0;
    for (int p = 0; p < world_size; ++p) {
      int rows = block_rows_per_proc[p] * kBlockSize;
      if (p == world_size - 1) {
        rows += height_remainder;
      }
      rows_per_proc[p] = rows;
      row_displs[p] = pixel_offset;
      pixel_offset += rows;
    }
  }

  auto compute_halo = [&](int p) -> std::pair<int, int> {
    if (rows_per_proc[p] == 0) {
      return {0, 0};
    }
    int top = (row_displs[p] > 0) ? 1 : 0;
    int bot = (row_displs[p] + rows_per_proc[p] < height) ? 1 : 0;
    return {top, bot};
  };

  const int local_rows = rows_per_proc[rank];
  const int local_block_rows = block_rows_per_proc[rank];
  auto [halo_top, halo_bottom] = compute_halo(rank);
  const int buffer_height = local_rows + halo_top + halo_bottom;

  std::vector<uint8_t> local_input(static_cast<size_t>(buffer_height) * width * 3);
  std::vector<int> scatter_counts(world_size);
  std::vector<int> scatter_displs(world_size);
  for (int p = 0; p < world_size; ++p) {
    auto [p_top, p_bot] = compute_halo(p);
    int p_buffer_rows = rows_per_proc[p] + p_top + p_bot;
    scatter_counts[p] = p_buffer_rows * width * 3;
    scatter_displs[p] = (row_displs[p] - p_top) * width * 3;
  }

  const uint8_t *send_buf = (rank == 0) ? std::get<2>(GetInput()).data() : nullptr;
  MPI_Scatterv(send_buf, scatter_counts.data(), scatter_displs.data(), MPI_UNSIGNED_CHAR,
               local_input.data(), static_cast<int>(local_input.size()), MPI_UNSIGNED_CHAR,
               0, MPI_COMM_WORLD);

  std::vector<uint8_t> local_output(static_cast<size_t>(local_rows) * width * 3);

  if (local_rows > 0) {
    int num_threads = ppc::util::GetNumThreads();
    if (num_threads <= 0) {
      num_threads = 1;
    }
    if (num_threads > local_rows) {
      num_threads = local_rows;
    }

    const bool is_last = (rank == world_size - 1);
    const int bottom_row_start = local_block_rows * kBlockSize;
    const int start_col_tail = num_col_blocks * kBlockSize;

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    auto processing = [&](int current_part) {
      int left_border_r = (local_block_rows * current_part) / num_threads;
      int right_border_r = (local_block_rows * (current_part + 1)) / num_threads;

      for (int bi = left_border_r; bi < right_border_r; ++bi) {
        for (int bj = 0; bj < num_col_blocks; ++bj) {
          ProcessFullBlock(local_input, local_output, width, buffer_height, halo_top,
                           bi * kBlockSize, bj * kBlockSize);
        }
      }

      if (width_has_remainder) {
        for (int bi = left_border_r; bi < right_border_r; ++bi) {
          ProcessPartBlock(local_input, local_output, width, local_rows, buffer_height, halo_top,
                           bi * kBlockSize, start_col_tail);
        }
      }

      if (is_last && height_remainder > 0) {
        int left_border_l = (num_col_blocks * current_part) / num_threads;
        int right_border_l = (num_col_blocks * (current_part + 1)) / num_threads;
        for (int bj = left_border_l; bj < right_border_l; ++bj) {
          ProcessPartBlock(local_input, local_output, width, local_rows, buffer_height, halo_top,
                           bottom_row_start, bj * kBlockSize);
        }
      }
    };

    for (int tid = 0; tid < num_threads; ++tid) {
      threads.emplace_back(processing, tid);
    }
    for (auto &th : threads) {
      th.join();
    }

    if (is_last && height_remainder > 0) {
      ProcessPartBlock(local_input, local_output, width, local_rows, buffer_height, halo_top,
                       bottom_row_start, start_col_tail);
    }
  }

  std::vector<int> recv_counts(world_size);
  std::vector<int> recv_displs(world_size);
  for (int p = 0; p < world_size; ++p) {
    recv_counts[p] = rows_per_proc[p] * width * 3;
    recv_displs[p] = row_displs[p] * width * 3;
  }

  std::vector<uint8_t> result(static_cast<size_t>(height) * width * 3);
  MPI_Gatherv(local_output.data(), static_cast<int>(local_output.size()), MPI_UNSIGNED_CHAR,
              result.data(), recv_counts.data(), recv_displs.data(), MPI_UNSIGNED_CHAR,
              0, MPI_COMM_WORLD);

  MPI_Bcast(result.data(), static_cast<int>(result.size()), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  GetOutput() = std::move(result);
  return true;
}

bool RomanovAGaussBlockALL::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_a_gauss_block