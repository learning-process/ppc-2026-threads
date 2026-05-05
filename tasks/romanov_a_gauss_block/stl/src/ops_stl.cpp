#include "romanov_a_gauss_block/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <tuple>
#include <vector>

#include "romanov_a_gauss_block/common/include/common.hpp"
#include "util/include/util.hpp"

namespace romanov_a_gauss_block {

RomanovAGaussBlockSTL::RomanovAGaussBlockSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool RomanovAGaussBlockSTL::ValidationImpl() {
  return std::get<0>(GetInput()) * std::get<1>(GetInput()) * 3 == static_cast<int>(std::get<2>(GetInput()).size());
}

bool RomanovAGaussBlockSTL::PreProcessingImpl() {
  return true;
}

namespace {

constexpr int kBlockSize = 32;

int ApplyKernel(const std::vector<uint8_t> &img, int row, int col, int channel, int width, int height,
                const std::array<std::array<int, 3>, 3> &kernel) {
  int sum = 0;
  for (size_t kr = 0; kr < 3; ++kr) {
    for (size_t kc = 0; kc < 3; ++kc) {
      int nr = row + static_cast<int>(kr) - 1;
      int nc = col + static_cast<int>(kc) - 1;
      if (nr >= 0 && nr < height && nc >= 0 && nc < width) {
        size_t idx = (((static_cast<size_t>(nr) * width) + nc) * 3) + channel;
        sum += (static_cast<int>(img[idx]) * kernel.at(kr).at(kc));
      }
    }
  }
  return sum;
}

void ProcessFullBlock(const std::vector<uint8_t> &initial_picture, std::vector<uint8_t> &result_picture, int width,
                      int height, int start_row, int start_col) {
  static constexpr std::array<std::array<int, 3>, 3> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

  for (int row = start_row; row < start_row + kBlockSize; ++row) {
    for (int col = start_col; col < start_col + kBlockSize; ++col) {
      for (int channel = 0; channel < 3; ++channel) {
        int sum = ApplyKernel(initial_picture, row, col, channel, width, height, kKernel);
        int result_value = (sum + 8) / 16;
        result_value = std::clamp(result_value, 0, 255);
        auto idx = ((static_cast<size_t>(row) * width + col) * 3) + channel;
        result_picture[idx] = static_cast<uint8_t>(result_value);
      }
    }
  }
}

void ProcessPartBlock(const std::vector<uint8_t> &initial_picture, std::vector<uint8_t> &result_picture, int width,
                      int height, int start_row, int start_col) {
  static constexpr std::array<std::array<int, 3>, 3> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

  const int end_row = std::min(height, start_row + kBlockSize);
  const int end_col = std::min(width, start_col + kBlockSize);

  for (int row = start_row; row < end_row; ++row) {
    for (int col = start_col; col < end_col; ++col) {
      for (int channel = 0; channel < 3; ++channel) {
        int sum = ApplyKernel(initial_picture, row, col, channel, width, height, kKernel);
        int result_value = (sum + 8) / 16;
        result_value = std::clamp(result_value, 0, 255);
        auto idx = ((static_cast<size_t>(row) * width + col) * 3) + channel;
        result_picture[idx] = static_cast<uint8_t>(result_value);
      }
    }
  }
}

}  // namespace

bool RomanovAGaussBlockSTL::RunImpl() {
  const int width = std::get<0>(GetInput());
  const int height = std::get<1>(GetInput());

  const std::vector<uint8_t> &initial_picture = std::get<2>(GetInput());
  std::vector<uint8_t> result_picture(static_cast<size_t>(height * width * 3));

  const int num_row_blocks = height / kBlockSize;
  const int num_col_blocks = width / kBlockSize;

  const int num_threads = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  auto processing = [&](int current_part) {
    int left_border_r = (num_row_blocks * current_part) / num_threads;
    int right_border_r = (num_row_blocks * (current_part + 1)) / num_threads;
    for (int bi = left_border_r; bi < right_border_r; ++bi) {
      for (int bj = 0; bj < num_col_blocks; ++bj) {
        ProcessFullBlock(initial_picture, result_picture, width, height, bi * kBlockSize, bj * kBlockSize);
      }
    }
    for (int bi = left_border_r; bi < right_border_r; ++bi) {
      ProcessPartBlock(initial_picture, result_picture, width, height, bi * kBlockSize, width - (width % kBlockSize));
    }
    int left_border_l = (num_col_blocks * current_part) / num_threads;
    int right_border_l = (num_col_blocks * (current_part + 1)) / num_threads;
    for (int bj = left_border_l; bj < right_border_l; ++bj) {
      ProcessPartBlock(initial_picture, result_picture, width, height, height - (height % kBlockSize), bj * kBlockSize);
    }
  };

  for (int tid = 0; tid < num_threads; tid++) {
    threads.emplace_back(processing, tid);
  }

  for (auto &th : threads) {
    th.join();
  }

  ProcessPartBlock(initial_picture, result_picture, width, height, height - (height % kBlockSize),
                   width - (width % kBlockSize));

  GetOutput() = result_picture;
  return true;
}

bool RomanovAGaussBlockSTL::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_a_gauss_block
