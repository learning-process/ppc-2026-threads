#include "chaschin_vladimir_linear_image_filtration_seq/seq/include/ops_seq.hpp"

#include <utility>
#include <vector>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"
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

  std::vector<float> temp(static_cast<std::vector<float>::size_type>(n * m));

  // ---------- горизонтальная фильтрация ----------
  for (int y = 0; y < m; ++y) {
    for (int x = 0; x < n; ++x) {
      float left = (x > 0) ? image[y * n + (x - 1)] : image[y * n + x];
      float center = image[y * n + x];
      float right = (x < n - 1) ? image[y * n + (x + 1)] : image[y * n + x];

      temp[y * n + x] = (left + 2.0f * center + right) / 4.0f;
    }
  }

  // ---------- вертикальная фильтрация ----------
  for (int x = 0; x < n; ++x) {
    for (int y = 0; y < m; ++y) {
      float top = (y > 0) ? temp[(y - 1) * n + x] : temp[y * n + x];
      float center = temp[y * n + x];
      float bottom = (y < m - 1) ? temp[(y + 1) * n + x] : temp[y * n + x];

      out[y * n + x] = (top + 2.0f * center + bottom) / 4.0f;
    }
  }

  return true;
}

bool ChaschinVLinearFiltrationSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace chaschin_v_linear_image_filtration_seq
