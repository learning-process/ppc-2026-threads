#include "kopilov_d_vertical_gauss_filter/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range2d.h>
#include <oneapi/tbb/parallel_for.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"

namespace kopilov_d_vertical_gauss_filter {

namespace {
const int kDivisor = 16;
const std::array<std::array<int, 3>, 3> kGaussKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t GetPixelMirroredTBB(const std::vector<uint8_t> &image, int x, int y, int width, int height) {
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
  auto idx = (static_cast<size_t>(new_y) * static_cast<size_t>(width)) + static_cast<size_t>(new_x);
  return image[idx];
}
}  // namespace

KopilovDVerticalGaussFilterTBB::KopilovDVerticalGaussFilterTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool KopilovDVerticalGaussFilterTBB::ValidationImpl() {
  const auto &in = GetInput();

  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  if (in.data.size() != static_cast<size_t>(in.width) * static_cast<size_t>(in.height)) {
    return false;
  }
  return true;
}

bool KopilovDVerticalGaussFilterTBB::PreProcessingImpl() {
  return true;
}

bool KopilovDVerticalGaussFilterTBB::RunImpl() {
  const auto &in = GetInput();
  int width = in.width;
  int height = in.height;
  const std::vector<uint8_t> &source_image = in.data;
  std::vector<uint8_t> destination_image(static_cast<size_t>(width) * static_cast<size_t>(height));

  tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&](const tbb::blocked_range2d<int> &range) {
    for (int j = range.rows().begin(); j != range.rows().end(); ++j) {
      for (int i = range.cols().begin(); i != range.cols().end(); ++i) {
        int pixel_sum = 0;

        for (size_t ky = 0; ky < 3; ++ky) {
          for (size_t kx = 0; kx < 3; ++kx) {
            int current_x = i + static_cast<int>(kx) - 1;
            int current_y = j + static_cast<int>(ky) - 1;
            pixel_sum +=
                kGaussKernel.at(ky).at(kx) * GetPixelMirroredTBB(source_image, current_x, current_y, width, height);
          }
        }

        auto out_idx = (static_cast<size_t>(j) * static_cast<size_t>(width)) + static_cast<size_t>(i);
        destination_image[out_idx] = static_cast<uint8_t>(pixel_sum / kDivisor);
      }
    }
  });

  GetOutput().width = width;
  GetOutput().height = height;
  GetOutput().data = std::move(destination_image);
  return true;
}

bool KopilovDVerticalGaussFilterTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kopilov_d_vertical_gauss_filter
