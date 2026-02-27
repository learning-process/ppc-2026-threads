#include "kazennova_a_sobel_operator/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace kazennova_a_sobel_operator {

const int kGx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};

const int kGy[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

SobelSeq::SobelSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SobelSeq::ValidationImpl() {
  const auto &in = GetInput();
  if (in.empty()) {
    return false;
  }
  size_t size = static_cast<size_t>(std::sqrt(in.size()));
  if (size * size != in.size()) {
    return false;
  }
  return true;
}

bool SobelSeq::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

static uint8_t GetPixel(const std::vector<uint8_t> &img, size_t size, int x, int y) {
  int size_int = static_cast<int>(size);

  if (x < 0) {
    x = 0;
  }
  if (x >= size_int) {
    x = size_int - 1;
  }
  if (y < 0) {
    y = 0;
  }
  if (y >= size_int) {
    y = size_int - 1;
  }

  size_t idx = static_cast<size_t>(y) * size + static_cast<size_t>(x);
  return img[idx];
}

bool SobelSeq::RunImpl() {
  const auto &in = GetInput();
  auto &out = GetOutput();

  size_t size = static_cast<size_t>(std::sqrt(in.size()));

  for (size_t y = 0; y < size; ++y) {
    for (size_t x = 0; x < size; ++x) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          uint8_t pixel = GetPixel(in, size, static_cast<int>(x) + kx, static_cast<int>(y) + ky);
          gx += static_cast<int>(pixel) * kGx[ky + 1][kx + 1];
          gy += static_cast<int>(pixel) * kGy[ky + 1][kx + 1];
        }
      }

      double magnitude_d = std::sqrt(static_cast<double>(gx * gx + gy * gy));
      int magnitude = static_cast<int>(std::round(magnitude_d));

      if (magnitude < 0) {
        magnitude = 0;
      }
      if (magnitude > 255) {
        magnitude = 255;
      }

      out[y * size + x] = static_cast<uint8_t>(magnitude);
    }
  }

  return true;
}

bool SobelSeq::PostProcessingImpl() {
  if (GetOutput().empty()) {
    return false;
  }
  return true;
}

}  // namespace kazennova_a_sobel_operator
