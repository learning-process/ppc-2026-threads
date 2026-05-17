#include "kolotukhin_a_gaussian_blur/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kolotukhin_a_gaussian_blur {

namespace {

void GetThreadRange(size_t tid, size_t total, size_t num_t, size_t &begin, size_t &end) {
  size_t chunk = total / num_t;
  begin = tid * chunk;
  end = (tid == num_t - 1) ? total : begin + chunk;
}

}  // namespace

KolotukhinAGaussinBlureSTL::KolotukhinAGaussinBlureSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KolotukhinAGaussinBlureSTL::ValidationImpl() {
  const auto &pixel_data = get<0>(GetInput());
  const auto img_width = get<1>(GetInput());
  const auto img_height = get<2>(GetInput());

  return static_cast<std::size_t>(img_height) * static_cast<std::size_t>(img_width) == pixel_data.size();
}

bool KolotukhinAGaussinBlureSTL::PreProcessingImpl() {
  const auto img_width = get<1>(GetInput());
  const auto img_height = get<2>(GetInput());

  GetOutput().assign(static_cast<std::size_t>(img_height) * static_cast<std::size_t>(img_width), 0);
  return true;
}

bool KolotukhinAGaussinBlureSTL::RunImpl() {
  const auto &pixel_data = get<0>(GetInput());
  const auto img_width = get<1>(GetInput());
  const auto img_height = get<2>(GetInput());

  const static std::array<std::array<int, 3>, 3> kKernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
  const static int kSum = 16;

  auto &output = GetOutput();
  const std::size_t total_pixels = static_cast<std::size_t>(img_height) * static_cast<std::size_t>(img_width);

  if (total_pixels == 0) {
    return false;
  }

  const int num_threads = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (int tid = 0; tid < num_threads; tid++) {
    threads.emplace_back([&, tid]() {
      std::size_t begin = 0;
      std::size_t end = 0;
      GetThreadRange(tid, total_pixels, num_threads, begin, end);

      for (std::size_t idx = begin; idx < end; idx++) {
        int row = static_cast<int>(idx / static_cast<std::size_t>(img_width));
        int col = static_cast<int>(idx % static_cast<std::size_t>(img_width));

        int acc = 0;
        for (int dy = -1; dy <= 1; dy++) {
          for (int dx = -1; dx <= 1; dx++) {
            std::uint8_t pixel = GetPixel(pixel_data, img_width, img_height, col + dx, row + dy);
            acc += kKernel.at(1 + dy).at(1 + dx) * static_cast<int>(pixel);
          }
        }
        output[idx] = static_cast<std::uint8_t>(acc / kSum);
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  return true;
}

std::uint8_t KolotukhinAGaussinBlureSTL::GetPixel(const std::vector<std::uint8_t> &pixel_data, int img_width,
                                                  int img_height, int pos_x, int pos_y) {
  std::size_t x = static_cast<std::size_t>(std::max(0, std::min(pos_x, img_width - 1)));
  std::size_t y = static_cast<std::size_t>(std::max(0, std::min(pos_y, img_height - 1)));
  return pixel_data[(y * static_cast<std::size_t>(img_width)) + x];
}

bool KolotukhinAGaussinBlureSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kolotukhin_a_gaussian_blur
