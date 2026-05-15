#include "iskhakov_d_vertical_gauss_filter/omp/include/ops_omp.hpp"

#include <omp.h>

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

uint8_t IskhakovDGetPixelMirrorOmp(const std::vector<uint8_t> &res, int col, int row, int width, int height) {
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
  return res[(row * width) + col];
}
}  // namespace

IskhakovDVerticalGaussFilterOMP::IskhakovDVerticalGaussFilterOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool IskhakovDVerticalGaussFilterOMP::ValidationImpl() {
  const auto &in = GetInput();

  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  if (in.data.size() != static_cast<size_t>(in.width) * static_cast<size_t>(in.height)) {
    return false;
  }
  return true;
}

bool IskhakovDVerticalGaussFilterOMP::PreProcessingImpl() {
  return true;
}

bool IskhakovDVerticalGaussFilterOMP::RunImpl() {
  const auto &in = GetInput();

  int width = in.width;
  int height = in.height;
  const std::vector<uint8_t> &matrix = in.data;
  std::vector<uint8_t> result(static_cast<size_t>(width) * static_cast<size_t>(height));

  int num_threads = ppc::util::GetNumThreads();
  omp_set_num_threads(num_threads);

#pragma omp parallel for default(none) shared(matrix, result, width, height, kGaussKernel, kDivConst)
  for (int horizontal_band = 0; horizontal_band < width; ++horizontal_band) {
    for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
      int sum = 0;

      sum += kGaussKernel[0][0] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band - 1, vertical_band - 1, width, height);
      sum += kGaussKernel[0][1] * IskhakovDGetPixelMirrorOmp(matrix, horizontal_band, vertical_band - 1, width, height);
      sum += kGaussKernel[0][2] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band + 1, vertical_band - 1, width, height);

      sum += kGaussKernel[1][0] * IskhakovDGetPixelMirrorOmp(matrix, horizontal_band - 1, vertical_band, width, height);
      sum += kGaussKernel[1][1] * IskhakovDGetPixelMirrorOmp(matrix, horizontal_band, vertical_band, width, height);
      sum += kGaussKernel[1][2] * IskhakovDGetPixelMirrorOmp(matrix, horizontal_band + 1, vertical_band, width, height);

      sum += kGaussKernel[2][0] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band - 1, vertical_band + 1, width, height);
      sum += kGaussKernel[2][1] * IskhakovDGetPixelMirrorOmp(matrix, horizontal_band, vertical_band + 1, width, height);
      sum += kGaussKernel[2][2] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band + 1, vertical_band + 1, width, height);

      result[(vertical_band * width) + horizontal_band] = static_cast<uint8_t>(sum / kDivConst);
    }
  }

  GetOutput().width = width;
  GetOutput().height = height;
  GetOutput().data = std::move(result);
  return true;
}

bool IskhakovDVerticalGaussFilterOMP::PostProcessingImpl() {
  return true;
}

}  // namespace iskhakov_d_vertical_gauss_filter
