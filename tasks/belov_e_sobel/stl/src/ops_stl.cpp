#include "belov_e_sobel/stl/include/ops_stl.hpp"

#include <thread>
#include <vector>
#include <cmath>
#include <algorithm>

#include "belov_e_sobel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace belov_e_sobel {

BelovESobelSTL::BelovESobelSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool BelovESobelSTL::ValidationImpl() {
  return !std::get<0>(GetInput()).empty() && (std::get<1>(GetInput()) > 0) && (std::get<2>(GetInput()) > 0);
}

bool BelovESobelSTL::PreProcessingImpl() {
  return true;
}

bool BelovESobelSTL::RunImpl() {
  const std::vector<uint8_t> &input = std::get<0>(GetInput());
  std::vector<uint8_t> &output = std::get<0>(GetOutput());
  int width = std::get<1>(GetInput());
  int height = std::get<2>(GetInput());

  auto get_px = [&](int x, int y) -> float {
    x = std::clamp(x, 0, width - 1);
    y = std::clamp(y, 0, height - 1);
    return static_cast<float>(input[(y * width) + x]);
  };

  unsigned int num_threads = ppc::util::GetNumThreads();

  std::vector<std::thread> threads;
  int rows_per_thread = height / num_threads;

  auto worker = [&](int start_y, int end_y) {
    for (int y = start_y; y < end_y; ++y) {
      for (int x = 0; x < width; ++x) {
        float gx = (-1 * get_px(x - 1, y - 1)) + (1 * get_px(x + 1, y - 1)) + (-2 * get_px(x - 1, y)) +
                   (2 * get_px(x + 1, y)) + (-1 * get_px(x - 1, y + 1)) + (1 * get_px(x + 1, y + 1));

        float gy = (-1 * get_px(x - 1, y - 1)) - (2 * get_px(x, y - 1)) - (1 * get_px(x + 1, y - 1)) +
                   (1 * get_px(x - 1, y + 1)) + (2 * get_px(x, y + 1)) + (1 * get_px(x + 1, y + 1));

        float magnitude = std::sqrt(gx * gx + gy * gy);
        output[y * width + x] = static_cast<uint8_t>(std::min(255.0f, magnitude));
      }
    }
  };

  for (unsigned int i = 0; i < num_threads; ++i) {
    int start_y = i * rows_per_thread;
    int end_y = (i == num_threads - 1) ? height : start_y + rows_per_thread;

    threads.emplace_back(worker, start_y, end_y);
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  return true;
}

bool BelovESobelSTL::PostProcessingImpl() {
  return !std::get<0>(GetOutput()).empty() && (std::get<1>(GetOutput()) > 0) && (std::get<2>(GetOutput()) > 0);
}

}  // namespace belov_e_sobel
