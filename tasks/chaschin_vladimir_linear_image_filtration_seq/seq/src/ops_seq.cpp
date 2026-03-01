#include "chaschin_vladimir_linear_image_filtration_seq/seq/include/ops_seq.hpp"

#include <numeric>
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
  int n = std::get<1>(in);
  int m = std::get<2>(in);

  auto &out = GetOutput();

  std::vector<float> temp(n * m);

  for (int y = 0; y < m; ++y) {
    temp[y * n + 0] = image[y * n + 0];
    temp[y * n + n - 1] = image[y * n + n - 1];

    for (int x = 1; x < n - 1; ++x) {
      temp[y * n + x] = (image[y * n + x - 1] + 2.0f * image[y * n + x] + image[y * n + x + 1]) * 0.25f;
    }
  }
  out.resize(n * m);

  for (int x = 0; x < n; ++x) {
    out[x] = temp[x];
    out[(m - 1) * n + x] = temp[(m - 1) * n + x];

    for (int y = 1; y < m - 1; ++y) {
      out[y * n + x] = (temp[(y - 1) * n + x] + 2.0f * temp[y * n + x] + temp[(y + 1) * n + x]) * 0.25f;
    }
  }
  return true;
}

bool ChaschinVLinearFiltrationSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace chaschin_v_linear_image_filtration_seq
