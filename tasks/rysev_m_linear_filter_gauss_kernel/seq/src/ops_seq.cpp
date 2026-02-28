#include "rysev_m_linear_filter_gauss_kernel/seq/include/ops_seq.hpp"

#include <stb/stb_image.h>

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "util/include/util.hpp"

namespace rysev_m_linear_filter_gauss_kernel {

RysevMGaussFilterSEQ::RysevMGaussFilterSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool RysevMGaussFilterSEQ::ValidationImpl() {
  return GetInput() >= 0 && GetOutput() == 0;
}

bool RysevMGaussFilterSEQ::PreProcessingImpl() {
  if (GetInput() == 0) {
    std::string abs_path =
        ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_rysev_m_linear_filter_gauss_kernel), "pic.ppm");
    int w, h, ch;
    unsigned char *data = stbi_load(abs_path.c_str(), &w, &h, &ch, STBI_rgb);
    if (data == nullptr) {
      return false;
    }
    width_ = w;
    height_ = h;
    channels_ = STBI_rgb;
    input_image_.assign(data, data + width_ * height_ * channels_);
    stbi_image_free(data);
  } else {
    int size = GetInput();
    width_ = height_ = size;
    channels_ = 3;
    input_image_.resize(width_ * height_ * channels_);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);
    for (auto &pixel : input_image_) {
      pixel = static_cast<uint8_t>(dist(gen));
    }
  }

  output_image_.resize(input_image_.size(), 0);
  return true;
}

bool RysevMGaussFilterSEQ::RunImpl() {
  const float kernel[3][3] = {
      {1.0f / 16, 2.0f / 16, 1.0f / 16}, {2.0f / 16, 4.0f / 16, 2.0f / 16}, {1.0f / 16, 2.0f / 16, 1.0f / 16}};

  int rows = height_;
  int cols = width_;
  int ch = channels_;

  for (int c = 0; c < ch; ++c) {
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        float sum = 0.0f;

        for (int di = -1; di <= 1; ++di) {
          for (int dj = -1; dj <= 1; ++dj) {
            int ni = i + di;
            int nj = j + dj;
            if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
              int idx = (ni * cols + nj) * ch + c;
              sum += static_cast<float>(input_image_[idx]) * kernel[di + 1][dj + 1];
            }
          }
        }

        int out_idx = (i * cols + j) * ch + c;
        output_image_[out_idx] = static_cast<uint8_t>(std::clamp(sum, 0.0f, 255.0f));
      }
    }
  }

  return true;
}

bool RysevMGaussFilterSEQ::PostProcessingImpl() {
  long long total = 0;
  for (uint8_t pixel : output_image_) {
    total += pixel;
  }
  GetOutput() = static_cast<int>(total);
  return true;
}

}  // namespace rysev_m_linear_filter_gauss_kernel
