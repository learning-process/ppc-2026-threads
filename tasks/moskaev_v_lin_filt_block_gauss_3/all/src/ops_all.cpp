#include "moskaev_v_lin_filt_block_gauss_3/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <thread>
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

void FilterPixelInBlock(const std::vector<uint8_t> &input_block, std::vector<uint8_t> &output_block, int block_w,
                        int channels, int row, int col, int ch) {
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

void FilterBlockRange(const std::vector<uint8_t> &input_block, std::vector<uint8_t> &output_block, int block_w,
                      int channels, int start_row, int end_row) {
  for (int row = start_row; row < end_row; ++row) {
    for (int col = 0; col < block_w; ++col) {
      for (int ch = 0; ch < channels; ++ch) {
        FilterPixelInBlock(input_block, output_block, block_w, channels, row, col, ch);
      }
    }
  }
}

void FilterBlock(const std::vector<uint8_t> &input_block, std::vector<uint8_t> &output_block, int block_w, int block_h,
                 int channels) {
  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  if (num_threads <= 1 || block_h < 2) {
    FilterBlockRange(input_block, output_block, block_w, channels, 0, block_h);
    return;
  }

  num_threads = std::min(num_threads, 8);
  num_threads = std::min(num_threads, block_h);
  int rows_per_thread = (block_h + num_threads - 1) / num_threads;
  std::vector<std::thread> threads;

  for (int tid = 0; tid < num_threads; ++tid) {
    int start = tid * rows_per_thread;
    int end = std::min(start + rows_per_thread, block_h);
    threads.emplace_back(FilterBlockRange, std::cref(input_block), std::ref(output_block), block_w, channels, start,
                         end);
  }

  for (auto &t : threads) {
    t.join();
  }
}

void ProcessOneBlock(int idx, int blocks_x, int width, int height, int channels, int block_size,
                     const std::vector<uint8_t> &image_data, std::vector<uint8_t> &output, int &output_offset) {
  int bx = idx % blocks_x;
  int by = idx / blocks_x;

  int block_x = bx * block_size;
  int block_y = by * block_size;
  int block_w = std::min(block_size, width - block_x);
  int block_h = std::min(block_size, height - block_y);
  int padded_w = block_w + 2;

  size_t input_size = static_cast<size_t>(padded_w) * static_cast<size_t>(block_h + 2) * static_cast<size_t>(channels);
  std::vector<uint8_t> input_block(input_size, 0);

  size_t output_size = static_cast<size_t>(block_w) * static_cast<size_t>(block_h) * static_cast<size_t>(channels);
  std::vector<uint8_t> output_block(output_size, 0);

  CopyBlockWithHalo(image_data, input_block, width, height, channels, block_x, block_y, block_w, block_h, padded_w);
  FilterBlock(input_block, output_block, block_w, block_h, channels);

  for (size_t i = 0; i < output_size; ++i) {
    output[output_offset + i] = output_block[i];
  }
  output_offset += static_cast<int>(output_size);
}

void BroadcastImageData(int rank, int &width, int &height, int &channels, std::vector<uint8_t> &image_data,
                        const InType &input) {
  if (rank == 0) {
    width = std::get<0>(input);
    height = std::get<1>(input);
    channels = std::get<2>(input);
    image_data = std::get<4>(input);
  }

  MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&channels, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int data_size = static_cast<int>(image_data.size());
  MPI_Bcast(&data_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    image_data.resize(data_size);
  }
  MPI_Bcast(image_data.data(), data_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
}

void ScatterBlocks(int rank, int num_procs, int total_blocks, std::vector<int> &local_blocks, int &local_cnt) {
  int per_proc = total_blocks / num_procs;
  int rem = total_blocks % num_procs;
  local_cnt = per_proc + (rank < rem ? 1 : 0);

  if (local_cnt <= 0) {
    local_blocks.clear();
    return;
  }

  std::vector<int> all(total_blocks);
  for (int i = 0; i < total_blocks; ++i) {
    all[i] = i;
  }

  std::vector<int> counts(num_procs);
  std::vector<int> displs(num_procs);
  int off = 0;
  for (int proc = 0; proc < num_procs; ++proc) {
    int cnt = per_proc + (proc < rem ? 1 : 0);
    counts[proc] = cnt;
    displs[proc] = off;
    off += cnt;
  }

  local_blocks.resize(local_cnt);
  MPI_Scatterv(all.data(), counts.data(), displs.data(), MPI_INT, local_blocks.data(), local_cnt, MPI_INT, 0,
               MPI_COMM_WORLD);
}

void ProcessBlockRange(const std::vector<int> &blocks, int start, int end, int blocks_x, int width, int height,
                       int channels, int block_size, const std::vector<uint8_t> &image_data,
                       std::vector<uint8_t> &output, int &output_offset) {
  for (int i = start; i < end; ++i) {
    ProcessOneBlock(blocks[i], blocks_x, width, height, channels, block_size, image_data, output, output_offset);
  }
}

void ProcessAssignedBlocksSequential(const std::vector<int> &local_blocks, int blocks_x, int width, int height,
                                     int channels, int block_size, const std::vector<uint8_t> &image_data,
                                     std::vector<uint8_t> &output) {
  int local_cnt = static_cast<int>(local_blocks.size());
  int total_bytes = 0;
  for (int i = 0; i < local_cnt; ++i) {
    int idx = local_blocks[i];
    int bx = idx % blocks_x;
    int by = idx / blocks_x;
    int block_x = bx * block_size;
    int block_y = by * block_size;
    int block_w = std::min(block_size, width - block_x);
    int block_h = std::min(block_size, height - block_y);
    total_bytes += block_w * block_h * channels;
  }
  output.resize(total_bytes);
  int output_offset = 0;
  ProcessBlockRange(local_blocks, 0, local_cnt, blocks_x, width, height, channels, block_size, image_data, output,
                    output_offset);
}

void ProcessBlocksInThread(int start, int blocks_in_thread, int blocks_x, int width, int height, int channels,
                           int block_size, const std::vector<uint8_t> &image_data, const std::vector<int> &local_blocks,
                           std::vector<uint8_t> &local_output) {
  int offset = 0;
  for (int i = start; i < start + blocks_in_thread; ++i) {
    int idx = local_blocks[i];
    int bx = idx % blocks_x;
    int by = idx / blocks_x;
    int block_x = bx * block_size;
    int block_y = by * block_size;
    int block_w = std::min(block_size, width - block_x);
    int block_h = std::min(block_size, height - block_y);
    int padded_w = block_w + 2;

    size_t input_size =
        static_cast<size_t>(padded_w) * static_cast<size_t>(block_h + 2) * static_cast<size_t>(channels);
    std::vector<uint8_t> input_block(input_size, 0);
    size_t output_size = static_cast<size_t>(block_w) * static_cast<size_t>(block_h) * static_cast<size_t>(channels);
    std::vector<uint8_t> output_block(output_size, 0);

    CopyBlockWithHalo(image_data, input_block, width, height, channels, block_x, block_y, block_w, block_h, padded_w);
    FilterBlock(input_block, output_block, block_w, block_h, channels);

    for (size_t j = 0; j < output_size; ++j) {
      local_output[offset + j] = output_block[j];
    }
    offset += static_cast<int>(output_size);
  }
}

void ProcessAssignedBlocksParallel(const std::vector<int> &local_blocks, int blocks_x, int width, int height,
                                   int channels, int block_size, const std::vector<uint8_t> &image_data,
                                   std::vector<uint8_t> &output) {
  int local_cnt = static_cast<int>(local_blocks.size());
  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  num_threads = std::min(num_threads, 8);
  num_threads = std::min(num_threads, local_cnt);
  int blocks_per_thread_base = local_cnt / num_threads;
  int blocks_remainder = local_cnt % num_threads;

  std::vector<std::vector<uint8_t>> thread_outputs(num_threads);
  std::vector<std::thread> threads;

  for (int tid = 0; tid < num_threads; ++tid) {
    int blocks_in_thread = blocks_per_thread_base + (tid < blocks_remainder ? 1 : 0);
    int start = (tid * blocks_per_thread_base) + std::min(tid, blocks_remainder);

    threads.emplace_back([&, tid, start, blocks_in_thread]() {
      int bytes_in_thread = 0;
      for (int i = start; i < start + blocks_in_thread; ++i) {
        int idx = local_blocks[i];
        int bx = idx % blocks_x;
        int by = idx / blocks_x;
        int block_x = bx * block_size;
        int block_y = by * block_size;
        int block_w = std::min(block_size, width - block_x);
        int block_h = std::min(block_size, height - block_y);
        bytes_in_thread += block_w * block_h * channels;
      }

      std::vector<uint8_t> local_output(bytes_in_thread);
      ProcessBlocksInThread(start, blocks_in_thread, blocks_x, width, height, channels, block_size, image_data,
                            local_blocks, local_output);
      thread_outputs[tid] = std::move(local_output);
    });
  }

  for (auto &t : threads) {
    t.join();
  }

  int total_bytes = 0;
  for (const auto &to : thread_outputs) {
    total_bytes += static_cast<int>(to.size());
  }
  output.resize(total_bytes);
  int pos = 0;
  for (const auto &to : thread_outputs) {
    std::ranges::copy(to, output.begin() + pos);
    pos += static_cast<int>(to.size());
  }
}

void ProcessAssignedBlocks(const std::vector<int> &local_blocks, int blocks_x, int width, int height, int channels,
                           int block_size, const std::vector<uint8_t> &image_data, std::vector<uint8_t> &output) {
  int local_cnt = static_cast<int>(local_blocks.size());
  if (local_cnt == 0) {
    output.clear();
    return;
  }

  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  if (num_threads <= 1 || local_cnt < 2) {
    ProcessAssignedBlocksSequential(local_blocks, blocks_x, width, height, channels, block_size, image_data, output);
  } else {
    ProcessAssignedBlocksParallel(local_blocks, blocks_x, width, height, channels, block_size, image_data, output);
  }
}

void GatherAndBroadcastResult(int rank, int num_procs, const std::vector<uint8_t> &output, OutType &out) {
  int send_count = static_cast<int>(output.size());

  std::vector<int> recv_counts(num_procs);
  MPI_Allgather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> displs(num_procs);
  int total_bytes = 0;
  for (int i = 0; i < num_procs; ++i) {
    displs[i] = total_bytes;
    total_bytes += recv_counts[i];
  }

  if (rank == 0) {
    out.resize(total_bytes);

    if (send_count > 0) {
      std::ranges::copy(output, out.begin());
    }

    for (int src = 1; src < num_procs; ++src) {
      if (recv_counts[src] > 0) {
        MPI_Recv(out.data() + displs[src], recv_counts[src], MPI_UNSIGNED_CHAR, src, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      }
    }
  } else {
    if (send_count > 0) {
      MPI_Send(output.data(), send_count, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
    }
  }

  int out_size = static_cast<int>(out.size());
  MPI_Bcast(&out_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    out.resize(out_size);
  }
  MPI_Bcast(out.data(), out_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
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
  if (rank_ != 0) {
    return true;
  }
  const auto &input = GetInput();
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

  BroadcastImageData(rank_, width, height, channels, image_data, GetInput());

  if (width == 0 || height == 0) {
    return false;
  }

  int blocks_x = (width + block_size_ - 1) / block_size_;
  int blocks_y = (height + block_size_ - 1) / block_size_;
  int total_blocks = blocks_x * blocks_y;

  if (total_blocks == 0) {
    return false;
  }

  std::vector<int> local_blocks;
  int local_cnt = 0;
  ScatterBlocks(rank_, num_procs_, total_blocks, local_blocks, local_cnt);

  std::vector<uint8_t> output;
  ProcessAssignedBlocks(local_blocks, blocks_x, width, height, channels, block_size_, image_data, output);

  GatherAndBroadcastResult(rank_, num_procs_, output, GetOutput());

  return true;
}

}  // namespace moskaev_v_lin_filt_block_gauss_3
