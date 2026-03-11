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

inline float HorizontalFilterAt(const std::vector<float> &img, int n, int x, int y) {
  const int idx = (y * n) + x;

  if (x == 0) {
    return ((2.F * img[idx]) + img[idx + 1]) / 3.F;
  }

  if (x == n - 1) {
    return (img[idx - 1] + (2.F * img[idx])) / 3.F;
  }

  return (img[idx - 1] + (2.F * img[idx]) + img[idx + 1]) / 4.F;
}
inline float VerticalFilterAt(const std::vector<float> &temp, int n, int m, int x, int y) {
  const int idx = (y * n) + x;

  if (y == 0) {
    return ((2.F * temp[idx]) + temp[idx + n]) / 3.F;
  }

  if (y == m - 1) {
    return (temp[idx - n] + (2.F * temp[idx])) / 3.F;
  }

  return (temp[idx - n] + (2.F * temp[idx]) + temp[idx + n]) / 4.F;
}

bool ChaschinVLinearFiltrationSEQ::RunImpl() {
  const auto &in = GetInput();
  const auto &image = std::get<0>(in);
  int n = std::get<1>(in);  // width
  int m = std::get<2>(in);  // height

  auto &out = GetOutput();
  out.resize(static_cast<std::vector<float>::size_type>(n) * static_cast<std::vector<float>::size_type>(m));

  std::vector<float> temp(static_cast<std::vector<float>::size_type>(n) *
                          static_cast<std::vector<float>::size_type>(m));

  // ---------- горизонтальная фильтрация ----------
  for (int yi = 0; yi < m; ++yi) {
    for (int xf = 0; xf < n; ++xf) {
      temp[(yi * n) + xf] = HorizontalFilterAt(image, n, xf, yi);
    }
  }

  for (int yi = 0; yi < m; ++yi) {
    for (int xy = 0; xy < n; ++xy) {
      out[(yi * n) + xy] = VerticalFilterAt(temp, n, m, xy, yi);
    }
  }

  return true;
}

bool ChaschinVLinearFiltrationSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace chaschin_v_linear_image_filtration_seq
