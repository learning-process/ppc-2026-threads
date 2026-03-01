#include "kopilov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"

namespace kopilov_d_vertical_gauss_filter {

namespace {
const int kDivConst = 16;
const std::array<std::array<int, 3>, 3> kGaussKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
const int kNumBands = 1;

uint8_t KopilovDGetPixelMirrorSeq(const std::vector<uint8_t> &res, int col, int row, int width, int height) {
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

KopilovDVerticalGaussFilterSEQ::KopilovDVerticalGaussFilterSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool KopilovDVerticalGaussFilterSEQ::ValidationImpl() {
  const auto &in = GetInput();

  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  if (in.data.size() != static_cast<size_t>(in.width) * static_cast<size_t>(in.height)) {
    return false;
  }
  return true;
}

bool KopilovDVerticalGaussFilterSEQ::PreProcessingImpl() {
  return true;
}

bool KopilovDVerticalGaussFilterSEQ::RunImpl() {
  const auto &in = GetInput();

  int width = in.width;
  int height = in.height;
  const std::vector<uint8_t> &matrix = in.data;
  std::vector<uint8_t> result(static_cast<size_t>(width) * static_cast<size_t>(height));

  int band_width = width / kNumBands;
  int remainder = width % kNumBands;
  int start_band = 0;

  for (int band = 0; band < kNumBands; ++band) {
    int cur_band_width = band_width + (band < remainder ? 1 : 0);
    int end_band = start_band + cur_band_width;

    for (int horizontal_band = start_band; horizontal_band < end_band; ++horizontal_band) {
      for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
        int sum = 0;

        sum += kGaussKernel[0][0] *
               KopilovDGetPixelMirrorSeq(matrix, horizontal_band - 1, vertical_band - 1, width, height);
        sum +=
            kGaussKernel[0][1] * KopilovDGetPixelMirrorSeq(matrix, horizontal_band, vertical_band - 1, width, height);
        sum += kGaussKernel[0][2] *
               KopilovDGetPixelMirrorSeq(matrix, horizontal_band + 1, vertical_band - 1, width, height);

        sum +=
            kGaussKernel[1][0] * KopilovDGetPixelMirrorSeq(matrix, horizontal_band - 1, vertical_band, width, height);
        sum += kGaussKernel[1][1] * KopilovDGetPixelMirrorSeq(matrix, horizontal_band, vertical_band, width, height);
        sum +=
            kGaussKernel[1][2] * KopilovDGetPixelMirrorSeq(matrix, horizontal_band + 1, vertical_band, width, height);

        sum += kGaussKernel[2][0] *
               KopilovDGetPixelMirrorSeq(matrix, horizontal_band - 1, vertical_band + 1, width, height);
        sum +=
            kGaussKernel[2][1] * KopilovDGetPixelMirrorSeq(matrix, horizontal_band, vertical_band + 1, width, height);
        sum += kGaussKernel[2][2] *
               KopilovDGetPixelMirrorSeq(matrix, horizontal_band + 1, vertical_band + 1, width, height);

        result[(vertical_band * width) + horizontal_band] = static_cast<uint8_t>(sum / kDivConst);
      }
    }
    start_band = end_band;
  }
  GetOutput().width = width;
  GetOutput().height = height;
  GetOutput().data = std::move(result);
  return true;
}

bool KopilovDVerticalGaussFilterSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kopilov_d_vertical_gauss_filter
