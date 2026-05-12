#include "nikolaev_d_block_linear_image_filtering/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "nikolaev_d_block_linear_image_filtering/common/include/common.hpp"

namespace nikolaev_d_block_linear_image_filtering {

NikolaevDBlockLinearImageFilteringTBB::NikolaevDBlockLinearImageFilteringTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool NikolaevDBlockLinearImageFilteringTBB::ValidationImpl() {
  const auto img_width = get<0>(GetInput());
  const auto img_height = get<1>(GetInput());
  const auto &pixel_data = get<2>(GetInput());

  return static_cast<std::size_t>(img_width) * static_cast<std::size_t>(img_height) * 3 == pixel_data.size();
}

bool NikolaevDBlockLinearImageFilteringTBB::PreProcessingImpl() {
  return true;
}

std::uint8_t NikolaevDBlockLinearImageFilteringTBB::GetPixel(const std::vector<uint8_t> &data, int w, int h, int nx,
                                                             int ny, int ch) {
  int ix = std::clamp(nx, 0, w - 1);
  int iy = std::clamp(ny, 0, h - 1);
  return data[((iy * w + ix) * 3) + ch];
}

std::uint8_t NikolaevDBlockLinearImageFilteringTBB::ApplyKernel(const std::vector<uint8_t> &src, int w, int h, int nx,
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

bool NikolaevDBlockLinearImageFilteringTBB::RunImpl() {
  const int width = std::get<0>(GetInput());
  const int height = std::get<1>(GetInput());
  const auto &src = std::get<2>(GetInput());

  auto &dst = GetOutput();
  dst.assign(src.size(), 0);

  tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&](const tbb::blocked_range2d<int> &r) {
    for (int ny = r.rows().begin(); ny < r.rows().end(); ++ny) {
      for (int nx = r.cols().begin(); nx < r.cols().end(); ++nx) {
        for (int ch = 0; ch < 3; ++ch) {
          dst[((ny * width + nx) * 3) + ch] = ApplyKernel(src, width, height, nx, ny, ch);
        }
      }
    }
  });

  return true;
}

bool NikolaevDBlockLinearImageFilteringTBB::PostProcessingImpl() {
  return true;
}

}  // namespace nikolaev_d_block_linear_image_filtering
