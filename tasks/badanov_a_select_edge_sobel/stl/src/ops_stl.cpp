#include "badanov_a_select_edge_sobel/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <thread>
#include <vector>

#include "badanov_a_select_edge_sobel/common/include/common.hpp"

namespace badanov_a_select_edge_sobel {

BadanovASelectEdgeSobelSTL::BadanovASelectEdgeSobelSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool BadanovASelectEdgeSobelSTL::ValidationImpl() {
  const auto &input = GetInput();
  return !input.empty();
}

bool BadanovASelectEdgeSobelSTL::PreProcessingImpl() {
  const auto &input = GetInput();

  width_ = static_cast<int>(std::sqrt(input.size()));
  height_ = width_;

  if (width_ * height_ != static_cast<int>(input.size())) {
    width_ = static_cast<int>(input.size());
    height_ = 1;
  }

  GetOutput() = std::vector<uint8_t>(input.size(), 0);

  return true;
}

void BadanovASelectEdgeSobelSTL::ApplySobelOperator(const std::vector<uint8_t> &input, std::vector<float> &magnitude,
                                                    float &max_magnitude) {
  const int height = height_;
  const int width = width_;

  if (height < 3 || width < 3) {
    max_magnitude = 0.0F;
    return;
  }

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  int rows_to_process = height - 2;
  unsigned int rows_per_thread = rows_to_process / num_threads;
  rows_per_thread = std::max<unsigned int>(rows_per_thread, 1);
  num_threads = (rows_to_process + rows_per_thread - 1) / rows_per_thread;

  std::vector<std::thread> threads;
  std::mutex max_mutex;
  float global_max = 0.0F;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    unsigned int start_row = 1 + (thread_idx * rows_per_thread);
    int end_row = static_cast<int>((thread_idx == num_threads - 1) ? (height - 1) : (start_row + rows_per_thread));

    threads.emplace_back([&, start_row, end_row]() {
      float local_max = 0.0F;

      for (int row = static_cast<int>(start_row); row < end_row; ++row) {
        for (int col = 1; col < width - 1; ++col) {
          float gradient_x = 0.0F;
          float gradient_y = 0.0F;

          ComputeGradientAtPixel(input, row, col, gradient_x, gradient_y);

          const float magnitude_value = std::sqrt((gradient_x * gradient_x) + (gradient_y * gradient_y));
          const size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col);
          magnitude[idx] = magnitude_value;

          local_max = std::max(magnitude_value, local_max);
        }
      }

      std::scoped_lock lock(max_mutex);
      global_max = std::max(local_max, global_max);
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  max_magnitude = global_max;
}

void BadanovASelectEdgeSobelSTL::ComputeGradientAtPixel(const std::vector<uint8_t> &input, int row, int col,
                                                        float &gradient_x, float &gradient_y) const {
  gradient_x = 0.0F;
  gradient_y = 0.0F;

  for (int kernel_row = -1; kernel_row <= 1; ++kernel_row) {
    for (int kernel_col = -1; kernel_col <= 1; ++kernel_col) {
      const size_t pixel_index =
          (static_cast<size_t>(row + kernel_row) * static_cast<size_t>(width_)) + static_cast<size_t>(col + kernel_col);
      const uint8_t pixel = input[pixel_index];

      const int kx_idx = kernel_row + 1;
      const int ky_idx = kernel_col + 1;
      const int kernel_x_value = kKernelX.at(static_cast<size_t>(kx_idx)).at(static_cast<size_t>(ky_idx));
      const int kernel_y_value = kKernelY.at(static_cast<size_t>(kx_idx)).at(static_cast<size_t>(ky_idx));

      gradient_x += static_cast<float>(pixel) * static_cast<float>(kernel_x_value);
      gradient_y += static_cast<float>(pixel) * static_cast<float>(kernel_y_value);
    }
  }
}

void BadanovASelectEdgeSobelSTL::ApplyThreshold(const std::vector<float> &magnitude, float max_magnitude,
                                                std::vector<uint8_t> &output) const {
  if (max_magnitude <= 0.0F) {
    std::ranges::fill(output, 0);
    return;
  }

  const float scale = 255.0F / max_magnitude;
  const size_t size = magnitude.size();
  const int threshold = threshold_;

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  size_t chunk_size = (size + num_threads - 1) / num_threads;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t start = thread_idx * chunk_size;
    size_t end = (thread_idx == num_threads - 1) ? size : (start + chunk_size);

    threads.emplace_back([&, start, end]() {
      for (size_t i = start; i < end; ++i) {
        output[i] = (magnitude[i] * scale > static_cast<float>(threshold)) ? 255 : 0;
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

bool BadanovASelectEdgeSobelSTL::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (height_ < 3 || width_ < 3) {
    output = input;
    return true;
  }

  std::vector<float> magnitude(input.size(), 0.0F);
  float max_magnitude = 0.0F;

  ApplySobelOperator(input, magnitude, max_magnitude);
  ApplyThreshold(magnitude, max_magnitude, output);

  return true;
}

bool BadanovASelectEdgeSobelSTL::PostProcessingImpl() {
  return true;
}

}  // namespace badanov_a_select_edge_sobel
