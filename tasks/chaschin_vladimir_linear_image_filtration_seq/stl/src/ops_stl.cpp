#include "chaschin_vladimir_linear_image_filtration_seq/stt/include/ops_stl.hpp"

#include <execution>
#include <utility>
#include <vector>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"

namespace chaschin_v_linear_image_filtration_stt {

ChaschinVLinearFiltrationSTL::ChaschinVLinearFiltrationSTL(const chaschin_v_linear_image_filtration_seq::InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  auto in_copy = in;
  GetInput() = std::move(in_copy);
  this->GetOutput().clear();
}

bool ChaschinVLinearFiltrationSTL::ValidationImpl() {
  const auto &in = GetInput();
  const auto &image = std::get<0>(in);
  return !image.empty();
}

bool ChaschinVLinearFiltrationSTL::PreProcessingImpl() {
  return true;
}

inline float HorizontalFilterAtSTL(const std::vector<float> &img, int n, int x, int y) {
  const int idx = (y * n) + x;
  if (x == 0) {
    return ((2.F * img[idx]) + img[idx + 1]) / 3.F;
  }
  if (x == n - 1) {
    return (img[idx - 1] + (2.F * img[idx])) / 3.F;
  }
  return (img[idx - 1] + (2.F * img[idx]) + img[idx + 1]) / 4.F;
}

inline float VerticalFilterAtSTL(const std::vector<float> &temp, int n, int m, int x, int y) {
  const int idx = (y * n) + x;
  if (y == 0) {
    return ((2.F * temp[idx]) + temp[idx + n]) / 3.F;
  }
  if (y == m - 1) {
    return (temp[idx - n] + (2.F * temp[idx])) / 3.F;
  }
  return (temp[idx - n] + (2.F * temp[idx]) + temp[idx + n]) / 4.F;
}

bool ChaschinVLinearFiltrationSTL::RunImpl() {
  const auto &in = GetInput();
  const auto &image = std::get<0>(in);
  int n = std::get<1>(in);
  int m = std::get<2>(in);

  auto &out = GetOutput();
  out.resize(static_cast<size_t>(n) * m);

  std::vector<float> temp(static_cast<size_t>(n) * m);

  std::vector<int> rows(m);
  std::iota(rows.begin(), rows.end(), 0);

  // ---------- горизонтальная фильтрация ----------
  std::for_each(std::execution::par, rows.begin(), rows.end(), [&](int yi) {
    for (int yi = r.begin(); yi < r.end(); ++yi) {
      for (int xf = 0; xf < n; ++xf) {
        temp[(yi * n) + xf] = HorizontalFilterAtSTL(image, n, xf, yi);
      }
    }
  });

  // ---------- вертикальная фильтрация ----------
  std::for_each(std::execution::par, rows.begin(), rows.end(), [&](int yi) {
    for (int yi = r.begin(); yi < r.end(); ++yi) {
      for (int xy = 0; xy < n; ++xy) {
        out[(yi * n) + xy] = VerticalFilterAtSTL(temp, n, m, xy, yi);
      }
    }
  });

  return true;
}

bool ChaschinVLinearFiltrationSTL::PostProcessingImpl() {
  return true;
}

}  // namespace chaschin_v_linear_image_filtration_stt
