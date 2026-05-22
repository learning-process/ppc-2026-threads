#include "lopatin_a_sobel_operator/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <thread>
#include <vector>

#include "lopatin_a_sobel_operator/common/include/common.hpp"
#include "util/include/util.hpp"

namespace lopatin_a_sobel_operator {

const std::array<std::array<int, 3>, 3> kSobelX = {std::array<int, 3>{-1, 0, 1}, std::array<int, 3>{-2, 0, 2},
                                                   std::array<int, 3>{-1, 0, 1}};

const std::array<std::array<int, 3>, 3> kSobelY = {std::array<int, 3>{-1, -2, -1}, std::array<int, 3>{0, 0, 0},
                                                   std::array<int, 3>{1, 2, 1}};

void LopatinASobelOperatorSTL::RunSobel(const Image &img, std::size_t start, std::size_t end,
                                        lopatin_a_sobel_operator::OutType &output) {
  const auto &input_data = img.pixels;

  for (std::size_t j = start; j < end; ++j) {  // processing only pixels with a full 3 x 3 neighborhood size
    for (std::size_t i = 1; i < img.width - 1; ++i) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          std::uint8_t pixel = input_data[((j + ky) * img.width) + (i + kx)];
          gx += pixel * kSobelX.at(ky + 1).at(kx + 1);
          gy += pixel * kSobelY.at(ky + 1).at(kx + 1);
        }
      }

      auto magnitude = static_cast<int>(round(std::sqrt((gx * gx) + (gy * gy))));
      output[(j * img.width) + i] = (magnitude > img.threshold) ? magnitude : 0;
    }
  }
}

LopatinASobelOperatorSTL::LopatinASobelOperatorSTL(const InType &in) : h_(in.height), w_(in.width) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool LopatinASobelOperatorSTL::ValidationImpl() {
  const auto &input = GetInput();
  return h_ * w_ == input.pixels.size();
}

bool LopatinASobelOperatorSTL::PreProcessingImpl() {
  GetOutput().resize(h_ * w_);
  return true;
}

bool LopatinASobelOperatorSTL::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  const auto rows = input.height - 2;
  const auto n_threads = std::min(static_cast<std::size_t>(ppc::util::GetNumThreads()), rows);

  std::vector<std::thread> threads;
  threads.reserve(n_threads);

  const std::size_t portion = rows / n_threads;
  const std::size_t tail = rows % n_threads;

  for (std::size_t i = 0; i < n_threads; ++i) {
    std::size_t start = 1 + (i * portion);
    std::size_t end = 1 + ((i + 1) * portion);
    threads.emplace_back(RunSobel, std::cref(input), start, end, std::ref(output));
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  if (tail != 0) {
    std::size_t start = input.height - 1 - tail;
    RunSobel(input, start, input.height - 1, output);
  }

  return true;
}

bool LopatinASobelOperatorSTL::PostProcessingImpl() {
  return true;
}

}  // namespace lopatin_a_sobel_operator
