#include "iskhakov_d_vertical_gauss_filter/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "iskhakov_d_vertical_gauss_filter/common/include/common.hpp"
#include "util/include/util.hpp"

namespace iskhakov_d_vertical_gauss_filter {

namespace {
const int kDivConst = 16;
const std::array<std::array<int, 3>, 3> kGaussKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t IskhakovDGetPixelMirrorAll(const std::vector<uint8_t> &src, int col, int row, int width, int height) {
  if (col < 0) {
    col = -col - 1;
  } else if (col >= width) {
    col = (2 * width) - col - 1;
  }
  if (row < 0) {
    row = -row - 1;
  } else if (row >= height) {
    row = (2 * height) - row - 1;
  }
  return src[(row * width) + col];
}

void ProcessLocalBlock(const std::vector<uint8_t> &matrix, std::vector<uint8_t> &local_result, int width, int height,
                       int start_col, int end_col, int local_cols) {
#pragma omp parallel for default(none) \
    shared(matrix, local_result, width, height, start_col, end_col, local_cols, kGaussKernel, kDivConst)
  for (int horizontal_band = start_col; horizontal_band < end_col; ++horizontal_band) {
    const int local_col_idx = horizontal_band - start_col;
    for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
      int sum = 0;

      sum += kGaussKernel[0][0] *
             IskhakovDGetPixelMirrorAll(matrix, horizontal_band - 1, vertical_band - 1, width, height);
      sum += kGaussKernel[0][1] * IskhakovDGetPixelMirrorAll(matrix, horizontal_band, vertical_band - 1, width, height);
      sum += kGaussKernel[0][2] *
             IskhakovDGetPixelMirrorAll(matrix, horizontal_band + 1, vertical_band - 1, width, height);

      sum += kGaussKernel[1][0] * IskhakovDGetPixelMirrorAll(matrix, horizontal_band - 1, vertical_band, width, height);
      sum += kGaussKernel[1][1] * IskhakovDGetPixelMirrorAll(matrix, horizontal_band, vertical_band, width, height);
      sum += kGaussKernel[1][2] * IskhakovDGetPixelMirrorAll(matrix, horizontal_band + 1, vertical_band, width, height);

      sum += kGaussKernel[2][0] *
             IskhakovDGetPixelMirrorAll(matrix, horizontal_band - 1, vertical_band + 1, width, height);
      sum += kGaussKernel[2][1] * IskhakovDGetPixelMirrorAll(matrix, horizontal_band, vertical_band + 1, width, height);
      sum += kGaussKernel[2][2] *
             IskhakovDGetPixelMirrorAll(matrix, horizontal_band + 1, vertical_band + 1, width, height);

      local_result[(vertical_band * local_cols) + local_col_idx] = static_cast<uint8_t>(sum / kDivConst);
    }
  }
}

void CopyOwnBand(const std::vector<uint8_t> &local_result, std::vector<uint8_t> &global_result, int width, int height,
                 int start_col, int end_col, int local_cols) {
  for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
    for (int horizontal_band = start_col; horizontal_band < end_col; ++horizontal_band) {
      const int local_col_idx = horizontal_band - start_col;
      global_result[(vertical_band * width) + horizontal_band] =
          local_result[(vertical_band * local_cols) + local_col_idx];
    }
  }
}

void ReceiveOtherBands(std::vector<uint8_t> &global_result, int width, int height, int size, int cols_per_proc,
                       int remainder) {
  for (int sender_rank = 1; sender_rank < size; ++sender_rank) {
    const int sender_start_col = (sender_rank * cols_per_proc) + std::min(sender_rank, remainder);
    const int sender_cols = cols_per_proc + (sender_rank < remainder ? 1 : 0);
    if (sender_cols == 0) {
      continue;
    }
    std::vector<uint8_t> recv_buf(static_cast<size_t>(sender_cols) * height);
    MPI_Recv(recv_buf.data(), static_cast<int>(recv_buf.size()), MPI_UNSIGNED_CHAR, sender_rank, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
      for (int col = 0; col < sender_cols; ++col) {
        global_result[(vertical_band * width) + sender_start_col + col] = recv_buf[(vertical_band * sender_cols) + col];
      }
    }
  }
}

void SendOwnBand(const std::vector<uint8_t> &local_result) {
  if (!local_result.empty()) {
    MPI_Send(local_result.data(), static_cast<int>(local_result.size()), MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
  }
}

void GatherAndBroadcast(const std::vector<uint8_t> &local_result, std::vector<uint8_t> &global_result, int width,
                        int height, int rank, int size, int cols_per_proc, int remainder, int start_col, int end_col,
                        int local_cols) {
  if (rank == 0) {
    CopyOwnBand(local_result, global_result, width, height, start_col, end_col, local_cols);
    ReceiveOtherBands(global_result, width, height, size, cols_per_proc, remainder);
  } else {
    SendOwnBand(local_result);
  }
  MPI_Bcast(global_result.data(), static_cast<int>(global_result.size()), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
}

}  // namespace

IskhakovDVerticalGaussFilterALL::IskhakovDVerticalGaussFilterALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool IskhakovDVerticalGaussFilterALL::ValidationImpl() {
  const auto &in = GetInput();
  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  if (in.data.size() != static_cast<size_t>(in.width) * static_cast<size_t>(in.height)) {
    return false;
  }
  return true;
}

bool IskhakovDVerticalGaussFilterALL::PreProcessingImpl() {
  return true;
}

bool IskhakovDVerticalGaussFilterALL::RunImpl() {
  const auto &in = GetInput();

  const int width = in.width;
  const int height = in.height;
  const std::vector<uint8_t> &matrix = in.data;

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int cols_per_proc = width / size;
  const int remainder = width % size;
  const int start_col = (rank * cols_per_proc) + std::min(rank, remainder);
  const int end_col = start_col + cols_per_proc + (rank < remainder ? 1 : 0);
  const int local_cols = end_col - start_col;

  std::vector<uint8_t> local_result(static_cast<size_t>(local_cols) * height);

  omp_set_num_threads(ppc::util::GetNumThreads());

  ProcessLocalBlock(matrix, local_result, width, height, start_col, end_col, local_cols);

  std::vector<uint8_t> global_result(static_cast<size_t>(width) * height);
  GatherAndBroadcast(local_result, global_result, width, height, rank, size, cols_per_proc, remainder, start_col,
                     end_col, local_cols);

  GetOutput().width = width;
  GetOutput().height = height;
  GetOutput().data = std::move(global_result);

  return true;
}

bool IskhakovDVerticalGaussFilterALL::PostProcessingImpl() {
  return true;
}

}  // namespace iskhakov_d_vertical_gauss_filter
