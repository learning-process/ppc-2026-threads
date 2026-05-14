#include "morozov_n_sobels_filter/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "morozov_n_sobels_filter/common/include/common.hpp"
#include "util/include/util.hpp"

namespace morozov_n_sobels_filter {

MorozovNSobelsFilterSTL::MorozovNSobelsFilterSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;

  result_image_.height = in.height;
  result_image_.width = in.width;
  result_image_.pixels.resize(result_image_.height * result_image_.width, 0);
}

bool MorozovNSobelsFilterSTL::ValidationImpl() {
  const Image &input = GetInput();
  return (input.height == result_image_.height) && (input.width == result_image_.width) &&
         (input.pixels.size() == result_image_.pixels.size());
}

bool MorozovNSobelsFilterSTL::PreProcessingImpl() {
  return true;
}

bool MorozovNSobelsFilterSTL::RunImpl() {
  const Image &input = GetInput();

  const int k_num_threads = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;
  threads.reserve((k_num_threads));

  size_t start_row = 1;
  size_t end_row = input.height - 1;
  size_t total_rows = end_row - start_row;

  size_t rows_per_thread = total_rows / k_num_threads;
  size_t remaining_rows = total_rows % k_num_threads;

  size_t current_start = start_row;

  for (int i = 0; i < k_num_threads; i++) {
    size_t num_rows = rows_per_thread + (std::cmp_less(i, remaining_rows) ? 1 : 0);

    if (num_rows > 0) {
      threads.emplace_back([this, &input, current_start, num_rows]() { this->Filter(input, current_start, num_rows); });
    }

    current_start += num_rows;
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  GetOutput() = result_image_;
  return true;
}

bool MorozovNSobelsFilterSTL::PostProcessingImpl() {
  return true;
}

void MorozovNSobelsFilterSTL::Filter(const Image &img, size_t start_row, size_t num_rows) {
  size_t end_row = start_row + num_rows;

  for (size_t id_y = start_row; id_y < end_row; id_y++) {
    for (size_t id_x = 1; id_x < img.width - 1; id_x++) {
      size_t pixel_id = (id_y * img.width) + id_x;
      result_image_.pixels[pixel_id] = CalculateNewPixelColor(img, id_x, id_y);
    }
  }
}

uint8_t MorozovNSobelsFilterSTL::CalculateNewPixelColor(const Image &img, size_t x, size_t y) {
  constexpr int kRadX = 1;
  constexpr int kRadY = 1;
  constexpr size_t kZero = 0;

  int grad_x = 0;
  int grad_y = 0;

  for (int row_offset = -kRadY; row_offset <= kRadY; row_offset++) {
    for (int col_offset = -kRadX; col_offset <= kRadX; col_offset++) {
      size_t id_x = std::clamp(x + col_offset, kZero, img.width - 1);
      size_t id_y = std::clamp(y + row_offset, kZero, img.height - 1);
      size_t pixel_id = (id_y * img.width) + id_x;

      grad_x += img.pixels[pixel_id] * kKernelX.at(row_offset + kRadY).at(col_offset + kRadX);
      grad_y += img.pixels[pixel_id] * kKernelY.at(row_offset + kRadY).at(col_offset + kRadX);
    }
  }

  int gradient = static_cast<int>(std::sqrt((grad_x * grad_x) + (grad_y * grad_y)));
  gradient = std::clamp(gradient, 0, 255);

  return static_cast<uint8_t>(gradient);
}
}  // namespace morozov_n_sobels_filter
