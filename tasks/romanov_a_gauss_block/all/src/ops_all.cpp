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

namespace {

constexpr int kBlockSize = 32;

struct Distribution {
  std::vector<int> rows_per_proc;
  std::vector<int> row_displs;
  int halo_top{0};
  int halo_bottom{0};
  int buffer_height{0};
};

int ApplyKernel(const std::vector<uint8_t> &img, int row, int col, int channel, int width, int buffer_height,
                int halo_top, const std::array<std::array<int, 3>, 3> &kernel) {
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

void ProcessFullBlock(const std::vector<uint8_t> &input, std::vector<uint8_t> &output, int width, int buffer_height,
                      int halo_top, int start_row, int start_col) {
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

void ProcessPartBlock(const std::vector<uint8_t> &input, std::vector<uint8_t> &output, int width, int local_rows,
                      int buffer_height, int halo_top, int start_row, int start_col) {
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

Distribution BuildDistribution(int rank, int world_size, int height) {
  const int total_block_rows = height / kBlockSize;
  const int height_remainder = height % kBlockSize;

  std::vector<int> block_rows_per_proc(world_size);
  const int base_blocks = total_block_rows / world_size;
  const int extra_blocks = total_block_rows % world_size;
  for (int proc = 0; proc < world_size; ++proc) {
    block_rows_per_proc[proc] = base_blocks + (proc < extra_blocks ? 1 : 0);
  }

  Distribution dist;
  dist.rows_per_proc.resize(world_size);
  dist.row_displs.resize(world_size);
  int pixel_offset = 0;
  for (int proc = 0; proc < world_size; ++proc) {
    int rows = block_rows_per_proc[proc] * kBlockSize;
    if (proc == world_size - 1) {
      rows += height_remainder;
    }
    dist.rows_per_proc[proc] = rows;
    dist.row_displs[proc] = pixel_offset;
    pixel_offset += rows;
  }

  // halo для текущего ранга
  if (dist.rows_per_proc[rank] > 0) {
    dist.halo_top = (dist.row_displs[rank] > 0) ? 1 : 0;
    dist.halo_bottom = (dist.row_displs[rank] + dist.rows_per_proc[rank] < height) ? 1 : 0;
  }
  dist.buffer_height = dist.rows_per_proc[rank] + dist.halo_top + dist.halo_bottom;

  return dist;
}

std::pair<int, int> HaloFor(int proc, const Distribution &dist, int height) {
  if (dist.rows_per_proc[proc] == 0) {
    return {0, 0};
  }
  int top = (dist.row_displs[proc] > 0) ? 1 : 0;
  int bot = (dist.row_displs[proc] + dist.rows_per_proc[proc] < height) ? 1 : 0;
  return {top, bot};
}

void ScatterWithHalo(int rank, int world_size, int width, int height, const Distribution &dist,
                     const uint8_t *full_image, std::vector<uint8_t> &local_input) {
  std::vector<int> scatter_counts(world_size);
  std::vector<int> scatter_displs(world_size);
  for (int proc = 0; proc < world_size; ++proc) {
    auto [proc_top, proc_bot] = HaloFor(proc, dist, height);
    int proc_buffer_rows = dist.rows_per_proc[proc] + proc_top + proc_bot;
    scatter_counts[proc] = proc_buffer_rows * width * 3;
    scatter_displs[proc] = (dist.row_displs[proc] - proc_top) * width * 3;
  }

  const uint8_t *send_buf = (rank == 0) ? full_image : nullptr;
  MPI_Scatterv(send_buf, scatter_counts.data(), scatter_displs.data(), MPI_UNSIGNED_CHAR, local_input.data(),
               static_cast<int>(local_input.size()), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
}

void ProcessThreadShare(int current_part, int num_threads, int local_block_rows, int num_col_blocks, int local_rows,
                        int height_remainder, bool is_last, bool width_has_remainder, int width, int buffer_height,
                        int halo_top, const std::vector<uint8_t> &local_input, std::vector<uint8_t> &local_output) {
  const int start_col_tail = num_col_blocks * kBlockSize;
  const int bottom_row_start = local_block_rows * kBlockSize;

  int left_border_r = (local_block_rows * current_part) / num_threads;
  int right_border_r = (local_block_rows * (current_part + 1)) / num_threads;
  for (int bi = left_border_r; bi < right_border_r; ++bi) {
    for (int bj = 0; bj < num_col_blocks; ++bj) {
      ProcessFullBlock(local_input, local_output, width, buffer_height, halo_top, bi * kBlockSize, bj * kBlockSize);
    }
  }

  if (width_has_remainder) {
    for (int bi = left_border_r; bi < right_border_r; ++bi) {
      ProcessPartBlock(local_input, local_output, width, local_rows, buffer_height, halo_top, bi * kBlockSize,
                       start_col_tail);
    }
  }

  if (is_last && height_remainder > 0) {
    int left_border_l = (num_col_blocks * current_part) / num_threads;
    int right_border_l = (num_col_blocks * (current_part + 1)) / num_threads;
    for (int bj = left_border_l; bj < right_border_l; ++bj) {
      ProcessPartBlock(local_input, local_output, width, local_rows, buffer_height, halo_top, bottom_row_start,
                       bj * kBlockSize);
    }
  }
}

void RunLocal(int rank, int world_size, int width, int height, const Distribution &dist,
              const std::vector<uint8_t> &local_input, std::vector<uint8_t> &local_output) {
  const int local_rows = dist.rows_per_proc[rank];
  if (local_rows == 0) {
    return;
  }

  const int total_block_rows = height / kBlockSize;
  const int height_remainder = height % kBlockSize;
  const int num_col_blocks = width / kBlockSize;
  const bool width_has_remainder = (width % kBlockSize) != 0;
  const int local_block_rows =
      (rank < total_block_rows % world_size) ? ((total_block_rows / world_size) + 1) : (total_block_rows / world_size);
  const bool is_last = (rank == world_size - 1);

  int num_threads = std::max(1, ppc::util::GetNumThreads());
  num_threads = std::min(num_threads, local_rows);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (int tid = 0; tid < num_threads; ++tid) {
    threads.emplace_back([&, tid]() {
      ProcessThreadShare(tid, num_threads, local_block_rows, num_col_blocks, local_rows, height_remainder, is_last,
                         width_has_remainder, width, dist.buffer_height, dist.halo_top, local_input, local_output);
    });
  }
  for (auto &th : threads) {
    th.join();
  }

  if (is_last && height_remainder > 0) {
    ProcessPartBlock(local_input, local_output, width, local_rows, dist.buffer_height, dist.halo_top,
                     local_block_rows * kBlockSize, num_col_blocks * kBlockSize);
  }
}

void GatherAndBroadcast(int world_size, int width, int height, const Distribution &dist,
                        const std::vector<uint8_t> &local_output, std::vector<uint8_t> &result) {
  std::vector<int> recv_counts(world_size);
  std::vector<int> recv_displs(world_size);
  for (int proc = 0; proc < world_size; ++proc) {
    recv_counts[proc] = dist.rows_per_proc[proc] * width * 3;
    recv_displs[proc] = dist.row_displs[proc] * width * 3;
  }
  result.assign(static_cast<size_t>(height) * width * 3, 0);
  MPI_Gatherv(local_output.data(), static_cast<int>(local_output.size()), MPI_UNSIGNED_CHAR, result.data(),
              recv_counts.data(), recv_displs.data(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(result.data(), static_cast<int>(result.size()), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
}

}  // namespace

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
  return std::get<0>(GetInput()) * std::get<1>(GetInput()) * 3 == static_cast<int>(std::get<2>(GetInput()).size());
}

bool RomanovAGaussBlockALL::PreProcessingImpl() {
  return true;
}

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

  const Distribution dist = BuildDistribution(rank, world_size, height);

  std::vector<uint8_t> local_input(static_cast<size_t>(dist.buffer_height) * width * 3);
  const uint8_t *full_image = (rank == 0) ? std::get<2>(GetInput()).data() : nullptr;
  ScatterWithHalo(rank, world_size, width, height, dist, full_image, local_input);

  std::vector<uint8_t> local_output(static_cast<size_t>(dist.rows_per_proc[rank]) * width * 3);
  RunLocal(rank, world_size, width, height, dist, local_input, local_output);

  std::vector<uint8_t> result;
  GatherAndBroadcast(world_size, width, height, dist, local_output, result);
  GetOutput() = std::move(result);
  return true;
}

bool RomanovAGaussBlockALL::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_a_gauss_block
