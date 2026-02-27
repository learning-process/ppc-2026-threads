#include "kazennova_a_sobel_operator/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "kazennova_a_sobel_operator/common/include/common.hpp"

namespace {

uint8_t GetPixel(const std::vector<uint8_t> &img, size_t size, int x, int y) {
  const int size_int = static_cast<int>(size);
  x = std::max(x, 0);
  x = std::min(x, size_int - 1);
  y = std::max(y, 0);
  y = std::min(y, size_int - 1);
  const size_t idx = (static_cast<size_t>(y) * size) + static_cast<size_t>(x);
  return img[idx];
}

}  // namespace

namespace kazennova_a_sobel_operator {

SobelSeq::SobelSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SobelSeq::ValidationImpl() {
  const auto &in = GetInput();
  if (in.empty()) {
    return false;
  }
  const auto size = static_cast<size_t>(std::sqrt(in.size()));
  return (size * size == in.size());
}

bool SobelSeq::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool SobelSeq::RunImpl() {
  const auto &in = GetInput();
  auto &out = GetOutput();

  const auto size = static_cast<size_t>(std::sqrt(in.size()));

  const std::array<std::array<int, 3>, 3> k_gx = {{{{-1, 0, 1}}, {{-2, 0, 2}}, {{-1, 0, 1}}}};

  const std::array<std::array<int, 3>, 3> k_gy = {{{{-1, -2, -1}}, {{0, 0, 0}}, {{1, 2, 1}}}};

  for (size_t row = 0; row < size; ++row) {
    for (size_t col = 0; col < size; ++col) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          const uint8_t pixel = GetPixel(in, size, static_cast<int>(col) + kx, static_cast<int>(row) + ky);

          const size_t ky_idx = static_cast<size_t>(ky + 1);
          const size_t kx_idx = static_cast<size_t>(kx + 1);

          gx += static_cast<int>(pixel) * k_gx.at(ky_idx).at(kx_idx);
          gy += static_cast<int>(pixel) * k_gy.at(ky_idx).at(kx_idx);
        }
      }

      const double magnitude_d = std::sqrt(static_cast<double>(gx * gx) + static_cast<double>(gy * gy));
      int magnitude = static_cast<int>(std::round(magnitude_d));
      magnitude = std::clamp(magnitude, 0, 255);

      out[(row * size) + col] = static_cast<uint8_t>(magnitude);
    }
  }

  return true;
}

bool SobelSeq::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace kazennova_a_sobel_operator
