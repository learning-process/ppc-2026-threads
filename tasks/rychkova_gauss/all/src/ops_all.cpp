#include "rychkova_gauss/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "rychkova_gauss/common/include/common.hpp"
#include "util/include/util.hpp"

namespace rychkova_gauss {

namespace {
int Mirror(int x, int xmin, int xmax) {
  if (x < xmin) {
    return 1;
  }
  if (x >= xmax) {
    return xmax - 1;
  }
  return x;
};

Pixel ComputePixel(const Image &image, std::size_t x, std::size_t y, std::size_t width, std::size_t height) {
  Pixel result = {.R = 0, .G = 0, .B = 0};
  for (int shift_x = -1; shift_x < 2; shift_x++) {
    for (int shift_y = -1; shift_y < 2; shift_y++) {
      int xn = Mirror(static_cast<int>(x) + shift_x, 0, static_cast<int>(width));
      int yn = Mirror(static_cast<int>(y) + shift_y, 0, static_cast<int>(height));
      auto current = image[yn][xn];
      result.R += static_cast<uint8_t>(static_cast<double>(current.R) * kKernel[shift_x + 1][shift_y + 1]);
      result.G += static_cast<uint8_t>(static_cast<double>(current.G) * kKernel[shift_x + 1][shift_y + 1]);
      result.B += static_cast<uint8_t>(static_cast<double>(current.B) * kKernel[shift_x + 1][shift_y + 1]);
    }
  }
  return result;
}
}  // namespace

RychkovaGaussALL::RychkovaGaussALL(const InType &in) : loutput_({}), goutput_({}) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool RychkovaGaussALL::ValidationImpl() {
  if (GetInput().empty()) {
    return false;
  }
  const auto len = GetInput()[0].size();
  return std::ranges::all_of(GetInput(), [len](const auto &row) { return row.size() == len; });
}

bool RychkovaGaussALL::PreProcessingImpl() {
  const auto &image = GetInput();
  const auto width = image[0].size();
  const auto height = image.size();
  loutput_.resize(width * height * 3, 0);
  goutput_.resize(width * height * 3, 0);
  return true;
}

bool RychkovaGaussALL::RunImpl() {
  const auto &image = GetInput();
  const auto width = image[0].size();
  const auto height = image.size();
  GetOutput() = Image(height, std::vector<Pixel>(width, Pixel(0, 0, 0)));

  int n = 0;
  int idx = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  MPI_Comm_rank(MPI_COMM_WORLD, &idx);

  size_t row_per_proc = height / n;
  size_t remainder_rows = height % n;
  size_t start = (idx * row_per_proc) + std::min(static_cast<size_t>(idx), remainder_rows);
  size_t end = ((idx + 1) * row_per_proc) + std::min(static_cast<size_t>(idx + 1), remainder_rows);

  size_t local_rows = end - start;

  std::vector<uint8_t> local_output(local_rows * width * 3);

  std::vector<int> counts(n);
  std::vector<int> displs(n);
  for (int i = 0; i < n; i++) {
    size_t rows_i = row_per_proc + (static_cast<size_t>(i) < remainder_rows ? 1 : 0);
    counts[i] = static_cast<int>(rows_i * width * 3);
    displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
  }

#pragma omp parallel for shared(local_output, image, width, height, start, local_rows) default(none) \
    num_threads(ppc::util::GetNumThreads())
  for (size_t local_j = 0; local_j < local_rows; local_j++) {
    size_t global_j = start + local_j;
    for (size_t i = 0; i < width; i++) {
      auto px = ComputePixel(image, i, global_j, width, height);
      size_t flat_idx = ((local_j * width) + i) * 3;
      local_output[flat_idx] = px.R;
      local_output[flat_idx + 1] = px.G;
      local_output[flat_idx + 2] = px.B;
    }
  }

  MPI_Gatherv(local_output.data(), static_cast<int>(local_output.size()), MPI_UINT8_T, goutput_.data(), counts.data(),
              displs.data(), MPI_UINT8_T, 0, MPI_COMM_WORLD);

  MPI_Bcast(goutput_.data(), static_cast<int>(goutput_.size()), MPI_UINT8_T, 0, MPI_COMM_WORLD);

  auto &output = GetOutput();

#pragma omp parallel for collapse(2) shared(height, width, output) default(none) num_threads(ppc::util::GetNumThreads())
  for (size_t j = 0; j < height; j++) {
    for (size_t i = 0; i < width; i++) {
      size_t flat_idx = ((j *

                          width) +
                         i) *
                        3;
      output[j][i] = Pixel(goutput_[flat_idx], goutput_[flat_idx + 1], goutput_[flat_idx + 2]);
    }
  }
  return true;
}

bool RychkovaGaussALL::PostProcessingImpl() {
  return true;
}

}  // namespace rychkova_gauss
