#include "kopilov_d_vertical_gauss_filter/omp/include/ops_omp.hpp"

#include <omp.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kopilov_d_vertical_gauss_filter {

namespace {
const int kDivisor = 16;
const std::array<std::array<int, 3>, 3> kGaussKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t GetPixelMirroredOmp(const std::vector<uint8_t> &image, int x, int y, int width, int height) {
  int new_x = x;
  int new_y = y;

  if (new_x < 0) {
    new_x = -new_x - 1;
  } else if (new_x >= width) {
    new_x = (2 * width) - new_x - 1;
  }
  if (new_y < 0) {
    new_y = -new_y - 1;
  } else if (new_y >= height) {
    new_y = (2 * height) - new_y - 1;
  }
  return image[(new_y * width) + new_x];
}
}  // namespace

KopilovDVerticalGaussFilterOMP::KopilovDVerticalGaussFilterOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool KopilovDVerticalGaussFilterOMP::ValidationImpl() {
  const auto &in = GetInput();

  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  if (in.data.size() != static_cast<size_t>(in.width) * static_cast<size_t>(in.height)) {
    return false;
  }
  return true;
}

bool KopilovDVerticalGaussFilterOMP::PreProcessingImpl() {
  return true;
}

bool KopilovDVerticalGaussFilterOMP::RunImpl() {
  const auto &in = GetInput();
  int width = in.width;
  int height = in.height;
  const std::vector<uint8_t> &source_image = in.data;
  std::vector<uint8_t> destination_image(static_cast<size_t>(width) * static_cast<size_t>(height));

  int num_threads = ppc::util::GetNumThreads();
  omp_set_num_threads(num_threads);

#pragma omp parallel for default(none) shared(source_image, destination_image, width, height, kGaussKernel, kDivisor)
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      int pixel_sum = 0;
      for (int kernel_y = -1; kernel_y <= 1; ++kernel_y) {
        for (int kernel_x = -1; kernel_x <= 1; ++kernel_x) {
          pixel_sum += kGaussKernel.at(kernel_y + 1).at(kernel_x + 1) *
                       GetPixelMirroredOmp(source_image, i + kernel_x, j + kernel_y, width, height);
        }
      }
      destination_image[(j * width) + i] = static_cast<uint8_t>(pixel_sum / kDivisor);
    }
  }

  GetOutput().width = width;
  GetOutput().height = height;
  GetOutput().data = std::move(destination_image);
  return true;
}

bool KopilovDVerticalGaussFilterOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kopilov_d_vertical_gauss_filter
