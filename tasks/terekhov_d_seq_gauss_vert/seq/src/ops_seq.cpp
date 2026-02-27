#include "terekhov_d_seq_gauss_vert/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "terekhov_d_seq_gauss_vert/common/include/common.hpp"

namespace terekhov_d_seq_gauss_vert {

TerekhovDGaussVertSEQ::TerekhovDGaussVertSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TerekhovDGaussVertSEQ::ValidationImpl() {
  const auto &input = GetInput();
  
  if (input.width <= 0 || input.height <= 0) {
    return false;
  }
  
  if (static_cast<int>(input.data.size()) != input.width * input.height) {
    return false;
  }
  
  return true;
}

bool TerekhovDGaussVertSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  width_ = input.width;
  height_ = input.height;
  
  GetOutput().width = width_;
  GetOutput().height = height_;
  GetOutput().data.resize(static_cast<size_t>(width_) * static_cast<size_t>(height_));
  
  int padded_width = width_ + 2;
  int padded_height = height_ + 2;
  padded_image_.resize(static_cast<size_t>(padded_width) * static_cast<size_t>(padded_height));
  
  for (int row = 0; row < padded_height; ++row) {
    for (int col = 0; col < padded_width; ++col) {
      int src_x = col - 1;
      int src_y = row - 1;
      
      if (src_x < 0) src_x = -src_x - 1;
      if (src_x >= width_) src_x = (2 * width_) - src_x - 1;
      if (src_y < 0) src_y = -src_y - 1;
      if (src_y >= height_) src_y = (2 * height_) - src_y - 1;
      
      size_t padded_idx = (static_cast<size_t>(row) * static_cast<size_t>(padded_width)) + static_cast<size_t>(col);
      size_t src_idx = (static_cast<size_t>(src_y) * static_cast<size_t>(width_)) + static_cast<size_t>(src_x);
      padded_image_[padded_idx] = input.data[src_idx];
    }
  }
  
  return true;
}

bool TerekhovDGaussVertSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  
  if (input.data.empty() || width_ <= 0 || height_ <= 0) {
    return false;
  }
  
  int padded_width = width_ + 2;
  
  const int num_bands = 4;
  int band_width = width_ / num_bands;
  if (band_width < 1) {
    band_width = 1;
  }
  
  for (int band = 0; band < num_bands; ++band) {
    int start_x = band * band_width;
    int end_x = (band == num_bands - 1) ? width_ : ((band + 1) * band_width);
    
    for (int row = 0; row < height_; ++row) {
      for (int col = start_x; col < end_x; ++col) {
        size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(width_)) + static_cast<size_t>(col);
        
        float sum = 0.0F;
        
        for (int ky = -1; ky <= 1; ++ky) {
          for (int kx = -1; kx <= 1; ++kx) {
            int px = col + kx + 1;
            int py = row + ky + 1;
            
            int kernel_idx = ((ky + 1) * 3) + (kx + 1);
            size_t padded_idx = (static_cast<size_t>(py) * static_cast<size_t>(padded_width)) + static_cast<size_t>(px);
            int pixel_value = padded_image_[padded_idx];
            
            sum += static_cast<float>(pixel_value) * kGaussKernel[static_cast<size_t>(kernel_idx)];
          }
        }
        
        output.data[idx] = static_cast<int>(std::lround(sum));
      }
    }
  }
  
  return true;
}

bool TerekhovDGaussVertSEQ::PostProcessingImpl() {
  return GetOutput().data.size() == (static_cast<size_t>(GetOutput().width) * static_cast<size_t>(GetOutput().height));
}

}  // namespace terekhov_d_seq_gauss_vert
