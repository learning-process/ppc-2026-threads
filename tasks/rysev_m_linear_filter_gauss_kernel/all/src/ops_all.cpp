#include "rysev_m_linear_filter_gauss_kernel/all/include/ops_all.hpp"

#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <thread>
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
    {KernelElement{.dr = -1, .dc = -1, .weight = 1.0F / 16},
     KernelElement{.dr = -1, .dc = 0, .weight = 2.0F / 16},
     KernelElement{.dr = -1, .dc = 1, .weight = 1.0F / 16},
     KernelElement{.dr = 0, .dc = -1, .weight = 2.0F / 16},
     KernelElement{.dr = 0, .dc = 0, .weight = 4.0F / 16},
     KernelElement{.dr = 0, .dc = 1, .weight = 2.0F / 16},
     KernelElement{.dr = 1, .dc = -1, .weight = 1.0F / 16},
     KernelElement{.dr = 1, .dc = 0, .weight = 2.0F / 16},
     KernelElement{.dr = 1, .dc = 1, .weight = 1.0F / 16}}};

// Вспомогательная функция для одного пикселя (та же, что в SEQ)
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

RysevMGaussFilterAll::RysevMGaussFilterAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool RysevMGaussFilterAll::ValidationImpl() {
  return GetInput() >= 0 && GetOutput() == 0;
}

bool RysevMGaussFilterAll::PreProcessingImpl() {
  if (GetInput() == 0) {
    std::string abs_path =
        ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_rysev_m_linear_filter_gauss_kernel), "pic.ppm");
    int w = 0, h = 0, ch = 0;
    unsigned char *data = stbi_load(abs_path.c_str(), &w, &h, &ch, STBI_rgb);
    if (data == nullptr) return false;
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
    for (auto &pixel : input_image_) pixel = static_cast<uint8_t>(dist(gen));
  }
  output_image_.resize(input_image_.size(), 0);
  return true;
}

void RysevMGaussFilterAll::ProcessBlockSequential(int start_row, int end_row) {
  for (int channel = 0; channel < channels_; ++channel) {
    for (int row = start_row; row < end_row; ++row) {
      for (int col = 0; col < width_; ++col) {
        float sum = ComputePixelValue(row, col, channel, height_, width_, channels_, input_image_);
        std::size_t out_idx = (static_cast<std::size_t>(row) * width_) + col;
        out_idx = (out_idx * channels_) + channel;
        output_image_[out_idx] = static_cast<uint8_t>(std::clamp(sum, 0.0F, 255.0F));
      }
    }
  }
}

bool RysevMGaussFilterAll::RunImpl() {
  int rows = height_;
  if (rows == 0) return true;

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) num_threads = 1;
  if (static_cast<int>(num_threads) > rows) num_threads = rows;

  std::vector<std::thread> workers;
  int rows_per_thread = rows / num_threads;
  int remainder = rows % num_threads;

  int start_row = 0;
  for (unsigned int t = 0; t < num_threads; ++t) {
    int end_row = start_row + rows_per_thread + (t < remainder ? 1 : 0);
    workers.emplace_back(&RysevMGaussFilterAll::ProcessBlockSequential, this, start_row, end_row);
    start_row = end_row;
  }

  for (auto &worker : workers) {
    if (worker.joinable()) worker.join();
  }

  return true;
}

bool RysevMGaussFilterAll::PostProcessingImpl() {
  int64_t total = 0;
  for (uint8_t pixel : output_image_) total += pixel;
  GetOutput() = static_cast<int>(total);
  return true;
}

}  // namespace rysev_m_linear_filter_gauss_kernel