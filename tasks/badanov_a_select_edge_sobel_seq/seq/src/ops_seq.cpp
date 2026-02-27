#include "badanov_a_select_edge_sobel_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>
#include <algorithm>

#include "badanov_a_select_edge_sobel_seq/common/include/common.hpp"
#include "util/include/util.hpp"

namespace badanov_a_select_edge_sobel_seq {

BadanovASelectEdgeSobelSEQ::BadanovASelectEdgeSobelSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool BadanovASelectEdgeSobelSEQ::ValidationImpl() {
  const auto& input = GetInput();
  return !input.empty();
}

bool BadanovASelectEdgeSobelSEQ::PreProcessingImpl() {
  const auto& input = GetInput();
  
  width_ = static_cast<int>(std::sqrt(input.size()));
  height_ = width_;
  
  if (width_ * height_ != static_cast<int>(input.size())) {
    width_ = static_cast<int>(input.size());
    height_ = 1;
  }
  
  GetOutput() = std::vector<uint8_t>(input.size(), 0);
  
  return true;
}

bool BadanovASelectEdgeSobelSEQ::RunImpl() {
  const auto& input = GetInput();
  auto& output = GetOutput();
  
  if (height_ < 3 || width_ < 3) {
    output = input;
    return true;
  }
  
  std::vector<float> magnitude(input.size(), 0.0f);
  float max_mag = 0.0f;
  
  for (int y = 1; y < height_ - 1; ++y) {
    for (int x = 1; x < width_ - 1; ++x) {
      float gx = 0.0f, gy = 0.0f;
      size_t idx = y * width_ + x;
      
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          uint8_t pixel = input[(y + ky) * width_ + (x + kx)];
          gx += static_cast<float>(pixel) * static_cast<float>(KERNEL_X[ky + 1][kx + 1]);
          gy += static_cast<float>(pixel) * static_cast<float>(KERNEL_Y[ky + 1][kx + 1]);
        }
      }
      
      float mag = std::sqrt(gx * gx + gy * gy);
      magnitude[idx] = mag;
      if (mag > max_mag) max_mag = mag;
    }
  }
  
  if (max_mag > 0) {
    float scale = 255.0f / max_mag;
    for (size_t i = 0; i < magnitude.size(); ++i) {
      output[i] = (magnitude[i] * scale > threshold_) ? 255 : 0;
    }
  } else {
    std::fill(output.begin(), output.end(), 0);
  }
  
  return true;
}

bool BadanovASelectEdgeSobelSEQ::PostProcessingImpl() {

  return true;
}

}  // namespace badanov_a_select_edge_sobel_seq