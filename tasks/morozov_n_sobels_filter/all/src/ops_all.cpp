#include "morozov_n_sobels_filter/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "morozov_n_sobels_filter/common/include/common.hpp"
#include "util/include/util.hpp"

namespace morozov_n_sobels_filter {

MorozovNSobelsFilterALL::MorozovNSobelsFilterALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;

  result_image_.height = in.height;
  result_image_.width = in.width;
  result_image_.pixels.resize(result_image_.height * result_image_.width, 0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);
}

bool MorozovNSobelsFilterALL::ValidationImpl() {
  const Image &input = GetInput();
  return (input.height == result_image_.height) && (input.width == result_image_.width) &&
         (input.pixels.size() == result_image_.pixels.size());
}

void MorozovNSobelsFilterALL::SplitRows(size_t global_rows, size_t proc_num, size_t &start, size_t &count) const {
  size_t base = global_rows / size_;
  size_t rem = global_rows % size_;

  count = base + (proc_num < rem ? 1 : 0);

  start = (proc_num * base) + std::min<size_t>(proc_num, rem);
}

void MorozovNSobelsFilterALL::SendImageDataFromZeroProc(const Image &global, size_t halo) {
  for (int proc = 0; proc < size_; proc++) {
    size_t proc_start{};
    size_t proc_rows{};

    SplitRows(global.height, proc, proc_start, proc_rows);

    size_t top = (proc_start == 0) ? 0 : halo;
    size_t bottom = (proc_start + proc_rows == global.height) ? 0 : halo;
    size_t offset = proc_start - top;
    size_t total = proc_rows + top + bottom;

    const auto *ptr = global.pixels.data() + (offset * global.width);

    if (proc == 0) {
      std::copy(ptr, ptr + (total * global.width), local_image_.pixels.begin());
    } else {
      MPI_Send(ptr, static_cast<int>(total * global.width), MPI_UNSIGNED_CHAR, proc, 0, MPI_COMM_WORLD);
    }
  }
}

void MorozovNSobelsFilterALL::CollectResult() {
  const int halo = 1;

  size_t top = (start_row_ == 0) ? 0 : halo;

  auto *send_ptr = result_image_.pixels.data() + (top * result_image_.width);

  size_t send_count = local_rows_ * result_image_.width;

  if (rank_ == 0) {
    Image global;

    global.width = result_image_.width;
    global.height = GetInput().height;
    global.pixels.resize(global.width * global.height);

    std::copy(send_ptr, send_ptr + send_count, global.pixels.begin());

    for (int proc = 1; proc < size_; proc++) {
      size_t start{};
      size_t rows{};

      SplitRows(global.height, proc, start, rows);

      MPI_Recv(global.pixels.data() + (start * global.width), static_cast<int>(rows * global.width), MPI_UNSIGNED_CHAR,
               proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    GetOutput() = global;

  } else {
    MPI_Send(send_ptr, static_cast<int>(send_count), MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
  }
}

bool MorozovNSobelsFilterALL::PreProcessingImpl() {
  const int halo = 1;

  const Image &global = GetInput();

  size_t width{};
  size_t height{};

  if (rank_ == 0) {
    width = global.width;
    height = global.height;
  }

  MPI_Bcast(&width, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&height, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  SplitRows(height, rank_, start_row_, local_rows_);

  size_t top_halo = (start_row_ == 0) ? 0 : halo;

  size_t bottom_halo = (start_row_ + local_rows_ == height) ? 0 : halo;

  size_t send_rows = local_rows_ + top_halo + bottom_halo;

  local_image_.width = width;
  local_image_.height = send_rows;
  local_image_.pixels.resize(send_rows * width);

  if (rank_ == 0) {
    SendImageDataFromZeroProc(global, static_cast<size_t>(halo));
  } else {
    MPI_Recv(local_image_.pixels.data(), static_cast<int>(send_rows * width), MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }

  return true;
}

bool MorozovNSobelsFilterALL::RunImpl() {
  result_image_.width = local_image_.width;
  result_image_.height = local_image_.height;
  result_image_.pixels.resize(local_image_.pixels.size());

  constexpr size_t kBegin = 1;
  size_t end = local_image_.height - 1;

  Filter(local_image_, result_image_, kBegin, end);
  CollectResult();

  return true;
}

void MorozovNSobelsFilterALL::Filter(const Image &img, Image &local_result, size_t start_row, size_t end_row) {
#pragma omp parallel for schedule(static) default(none) shared(img, local_result, start_row, end_row) \
    num_threads(ppc::util::GetNumThreads())
  for (size_t id_y = start_row; id_y < end_row; id_y++) {
    for (size_t id_x = 1; id_x < img.width - 1; id_x++) {
      size_t pixel_id = (id_y * img.width) + id_x;
      local_result.pixels[pixel_id] = CalculateNewPixelColor(img, id_x, id_y);
    }
  }
}

uint8_t MorozovNSobelsFilterALL::CalculateNewPixelColor(const Image &img, size_t x, size_t y) {
  constexpr int kRadX = 1;
  constexpr int kRadY = 1;
  constexpr size_t kZero = 0;

  int grad_x = 0;
  int grad_y = 0;

  for (int row_offset = -kRadY; row_offset <= kRadY; row_offset++) {
    for (int col_offset = -kRadX; col_offset <= kRadX; col_offset++) {
      size_t id_x = std::clamp(x + col_offset, kZero, img.width - 1);
      size_t id_y = std::clamp(y + row_offset, kZero, img.height - 1);
      size_t pixel_id = (id_y * img.width) + id_x;

      grad_x += img.pixels[pixel_id] * kKernelX.at(row_offset + kRadY).at(col_offset + kRadX);
      grad_y += img.pixels[pixel_id] * kKernelY.at(row_offset + kRadY).at(col_offset + kRadX);
    }
  }

  int gradient = static_cast<int>(std::sqrt((grad_x * grad_x) + (grad_y * grad_y)));
  gradient = std::clamp(gradient, 0, 255);

  return static_cast<uint8_t>(gradient);
}

bool MorozovNSobelsFilterALL::PostProcessingImpl() {
  Image &out = GetOutput();

  if (rank_ != 0) {
    out.width = result_image_.width;

    out.height = GetInput().height;

    out.pixels.resize(out.width * out.height);
  }

  MPI_Bcast(out.pixels.data(), static_cast<int>(out.pixels.size()), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  return true;
}

}  // namespace morozov_n_sobels_filter
