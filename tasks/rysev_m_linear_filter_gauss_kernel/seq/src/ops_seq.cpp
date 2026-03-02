#include "rysev_m_linear_filter_gauss_kernel/seq/include/ops_seq.hpp"

#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "rysev_m_linear_filter_gauss_kernel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace rysev_m_linear_filter_gauss_kernel {

namespace {
struct KernelElement {
  int dr;
  int dc;
  float weight;
};

const std::array<KernelElement, 9> kKernelElements = {
    {KernelElement{.dr = -1, .dc = -1, .weight = 1.0F / 16}, KernelElement{.dr = -1, .dc = 0, .weight = 2.0F / 16},
     KernelElement{.dr = -1, .dc = 1, .weight = 1.0F / 16}, KernelElement{.dr = 0, .dc = -1, .weight = 2.0F / 16},
     KernelElement{.dr = 0, .dc = 0, .weight = 4.0F / 16}, KernelElement{.dr = 0, .dc = 1, .weight = 2.0F / 16},
     KernelElement{.dr = 1, .dc = -1, .weight = 1.0F / 16}, KernelElement{.dr = 1, .dc = 0, .weight = 2.0F / 16},
     KernelElement{.dr = 1, .dc = 1, .weight = 1.0F / 16}}};

float ComputePixelValue(int row, int col, int channel, int rows, int cols, int channels,
                        const std::vector<uint8_t> &input) {
  float sum = 0.0F;
  for (const auto &ke : kKernelElements) {
    int nr = row + ke.dr;
    int nc = col + ke.dc;
    if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
      std::size_t idx = (static_cast<std::size_t>(nr) * cols) + nc;
      idx = (idx * channels) + channel;
      sum += static_cast<float>(input[idx]) * ke.weight;
    }
  }
  return sum;
}

}  // namespace

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
    int w = 0;
    int h = 0;
    int ch = 0;
    unsigned char *data = stbi_load(abs_path.c_str(), &w, &h, &ch, STBI_rgb);
    if (data == nullptr) {
      return false;
    }
    width_ = w;
    height_ = h;
    channels_ = STBI_rgb;
    std::ptrdiff_t total_pixels = static_cast<std::ptrdiff_t>(width_) * height_ * channels_;
    input_image_.assign(data, data + total_pixels);
    stbi_image_free(data);
  } else {
    int size = GetInput();
    width_ = height_ = size;
    channels_ = 3;
    std::size_t total_pixels = static_cast<std::size_t>(width_) * height_ * channels_;
    input_image_.resize(total_pixels);

    std::mt19937 gen(static_cast<unsigned int>(GetInput()));
    std::uniform_int_distribution<int> dist(0, 255);
    for (auto &pixel : input_image_) {
      pixel = static_cast<uint8_t>(dist(gen));
    }
  }

  output_image_.resize(input_image_.size(), 0);
  return true;
}

void RysevMGaussFilterSEQ::ApplyKernelToChannel(int channel, int rows, int cols) {
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      float sum = ComputePixelValue(row, col, channel, rows, cols, channels_, input_image_);
      std::size_t out_idx = (static_cast<std::size_t>(row) * cols) + col;
      out_idx = (out_idx * channels_) + channel;
      output_image_[out_idx] = static_cast<uint8_t>(std::clamp(sum, 0.0F, 255.0F));
    }
  }
}

bool RysevMGaussFilterSEQ::RunImpl() {
  int rows = height_;
  int cols = width_;
  int ch = channels_;

  for (int channel = 0; channel < ch; ++channel) {
    ApplyKernelToChannel(channel, rows, cols);
  }

  return true;
}

bool RysevMGaussFilterSEQ::PostProcessingImpl() {
  int64_t total = 0;
  for (uint8_t pixel : output_image_) {
    total += pixel;
  }
  GetOutput() = static_cast<int>(total);
  return true;
}

}  // namespace rysev_m_linear_filter_gauss_kernel
