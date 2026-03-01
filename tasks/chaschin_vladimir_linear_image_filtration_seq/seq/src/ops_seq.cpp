#include "chaschin_vladimir_linear_image_filtration_seq/seq/include/ops_seq.hpp"

#include <utility>
#include <vector>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"
#include "util/include/util.hpp"

namespace chaschin_v_linear_image_filtration_seq {

ChaschinVLinearFiltrationSEQ::ChaschinVLinearFiltrationSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  auto in_copy = in;
  GetInput() = std::move(in_copy);
  this->GetOutput().clear();
}

bool ChaschinVLinearFiltrationSEQ::ValidationImpl() {
  const auto &in = GetInput();

  const auto &image = std::get<0>(in);

  return !image.empty();
}

bool ChaschinVLinearFiltrationSEQ::PreProcessingImpl() {
  return true;
}

bool ChaschinVLinearFiltrationSEQ::RunImpl() {
  const auto &in = GetInput();
  const auto &image = std::get<0>(in);
  int n = std::get<1>(in);  // width
  int m = std::get<2>(in);  // height

  auto &out = GetOutput();
  out.resize(static_cast<std::vector<float>::size_type>(n * m));

  const float k[3][3] = {
      {1.f / 16.f, 2.f / 16.f, 1.f / 16.f}, {2.f / 16.f, 4.f / 16.f, 2.f / 16.f}, {1.f / 16.f, 2.f / 16.f, 1.f / 16.f}};

  for (int y = 0; y < m; ++y) {
    for (int x = 0; x < n; ++x) {
      float acc = 0.0f;
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          int ny = y + dy;
          int nx = x + dx;
          if (ny >= 0 && ny < m && nx >= 0 && nx < n) {
            acc += image[ny * n + nx] * k[dy + 1][dx + 1];
          }
        }
      }
      out[y * n + x] = acc;
    }
  }

  return true;
}

bool ChaschinVLinearFiltrationSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace chaschin_v_linear_image_filtration_seq
