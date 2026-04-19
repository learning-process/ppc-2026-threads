#include "lopatin_a_sobel_operator/seq/include/ops_seq.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "lopatin_a_sobel_operator/common/include/common.hpp"

namespace lopatin_a_sobel_operator {

const std::array<std::array<int, 3>, 3> kSobelX = {std::array<int, 3>{-1, 0, 1}, std::array<int, 3>{-2, 0, 2},
                                                   std::array<int, 3>{-1, 0, 1}};

const std::array<std::array<int, 3>, 3> kSobelY = {std::array<int, 3>{-1, -2, -1}, std::array<int, 3>{0, 0, 0},
                                                   std::array<int, 3>{1, 2, 1}};

LopatinASobelOperatorSEQ::LopatinASobelOperatorSEQ(const InType &in) : h_(in.height), w_(in.width) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool LopatinASobelOperatorSEQ::ValidationImpl() {
  const auto &input = GetInput();
  return h_ * w_ == input.pixels.size();
}

bool LopatinASobelOperatorSEQ::PreProcessingImpl() {
  GetOutput().resize(h_ * w_);
  return true;
}

bool LopatinASobelOperatorSEQ::RunImpl() {
  const auto &input = GetInput();
  const auto &input_data = input.pixels;
  auto &output = GetOutput();

  for (std::size_t j = 1; j < h_ - 1; ++j) {  // processing only pixels with a full 3 x 3 neighborhood size
    for (std::size_t i = 1; i < w_ - 1; ++i) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          std::uint8_t pixel = input_data[((j + ky) * w_) + (i + kx)];
          gx += pixel * kSobelX.at(ky + 1).at(kx + 1);
          gy += pixel * kSobelY.at(ky + 1).at(kx + 1);
        }
      }

      auto magnitude = static_cast<int>(round(std::sqrt((gx * gx) + (gy * gy))));
      output[(j * w_) + i] = (magnitude > input.threshold) ? magnitude : 0;
    }
  }

  return true;
}

bool LopatinASobelOperatorSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace lopatin_a_sobel_operator
