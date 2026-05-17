#include "nikolaev_d_block_linear_image_filtering/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <numeric>
#include <vector>

#include "nikolaev_d_block_linear_image_filtering/common/include/common.hpp"

namespace nikolaev_d_block_linear_image_filtering {

NikolaevDBlockLinearImageFilteringSTL::NikolaevDBlockLinearImageFilteringSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool NikolaevDBlockLinearImageFilteringSTL::ValidationImpl() {
  const auto img_width = std::get<0>(GetInput());
  const auto img_height = std::get<1>(GetInput());
  const auto &pixel_data = std::get<2>(GetInput());

  return static_cast<std::size_t>(img_width) * static_cast<std::size_t>(img_height) * 3 == pixel_data.size();
}

bool NikolaevDBlockLinearImageFilteringSTL::PreProcessingImpl() {
  return true;
}

std::uint8_t NikolaevDBlockLinearImageFilteringSTL::GetPixel(const std::vector<uint8_t> &data, int w, int h, int nx,
                                                             int ny, int ch) {
  int ix = std::clamp(nx, 0, w - 1);
  int iy = std::clamp(ny, 0, h - 1);
  return data[((iy * w + ix) * 3) + ch];
}

std::uint8_t NikolaevDBlockLinearImageFilteringSTL::ApplyKernel(const std::vector<uint8_t> &src, int w, int h, int nx,
                                                                int ny, int ch) {
  static constexpr std::array<std::array<int, 3>, 3> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
  int acc = 0;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      acc += GetPixel(src, w, h, nx + kx, ny + ky, ch) * kKernel.at(ky + 1).at(kx + 1);
    }
  }
  return static_cast<uint8_t>(std::clamp((acc + 8) / 16, 0, 255));
}

bool NikolaevDBlockLinearImageFilteringSTL::RunImpl() {
  const int width = std::get<0>(GetInput());
  const int height = std::get<1>(GetInput());
  const auto &src = std::get<2>(GetInput());

  auto &dst = GetOutput();
  dst.assign(src.size(), 0);

  std::vector<int> row_indices(height);
  std::iota(row_indices.begin(), row_indices.end(), 0);

  std::for_each(std::execution::par_unseq, row_indices.begin(), row_indices.end(), [&](int ny) {
    for (int nx = 0; nx < width; ++nx) {
      for (int ch = 0; ch < 3; ++ch) {
        dst[((ny * width + nx) * 3) + ch] = ApplyKernel(src, width, height, nx, ny, ch);
      }
    }
  });

  return true;
}

bool NikolaevDBlockLinearImageFilteringSTL::PostProcessingImpl() {
  return true;
}

}  // namespace nikolaev_d_block_linear_image_filtering
