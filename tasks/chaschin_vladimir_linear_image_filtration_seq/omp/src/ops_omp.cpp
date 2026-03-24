#include "chaschin_vladimir_linear_image_filtration_seq/omp/include/ops_omp.hpp"

#include <omp.h>

#include <utility>
#include <vector>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"
namespace chaschin_v_linear_image_filtration_omp {

ChaschinVLinearFiltrationOMP::ChaschinVLinearFiltrationOMP(const chaschin_v_linear_image_filtration::InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  auto in_copy = in;
  GetInput() = std::move(in_copy);
  this->GetOutput().clear();
}

bool ChaschinVLinearFiltrationOMP::ValidationImpl() {
  const auto &in = GetInput();
  const auto &image = std::get<0>(in);
  return !image.empty();
}

bool ChaschinVLinearFiltrationOMP::PreProcessingImpl() {
  return true;
}

inline float HorizontalFilterAtOMP(const std::vector<float> &img, int n, int x, int y) {
  const int idx = (y * n) + x;

  if (x == 0) {
    return ((2.F * img[idx]) + img[idx + 1]) / 3.F;
  }

  if (x == n - 1) {
    return (img[idx - 1] + (2.F * img[idx])) / 3.F;
  }

  return (img[idx - 1] + (2.F * img[idx]) + img[idx + 1]) / 4.F;
}

inline float VerticalFilterAtOMP(const std::vector<float> &temp, int n, int m, int x, int y) {
  const int idx = (y * n) + x;

  if (y == 0) {
    return ((2.F * temp[idx]) + temp[idx + n]) / 3.F;
  }

  if (y == m - 1) {
    return (temp[idx - n] + (2.F * temp[idx])) / 3.F;
  }

  return (temp[idx - n] + (2.F * temp[idx]) + temp[idx + n]) / 4.F;
}

bool ChaschinVLinearFiltrationOMP::RunImpl() {
  const auto &in = GetInput();
  const auto &image = std::get<0>(in);
  int n = std::get<1>(in);
  int m = std::get<2>(in);

  auto &out = GetOutput();
  out.resize(static_cast<std::vector<float>::size_type>(n) * static_cast<std::vector<float>::size_type>(m));

  std::vector<float> temp(static_cast<std::vector<float>::size_type>(n) *
                          static_cast<std::vector<float>::size_type>(m));

  // ---------- горизонтальная фильтрация ----------
#pragma omp parallel for default(none) shared(n, m, image, temp)
  for (int yi = 0; yi < m; ++yi) {
    for (int xf = 0; xf < n; ++xf) {
      temp[(yi * n) + xf] = HorizontalFilterAtOMP(image, n, xf, yi);
    }
  }

  // ---------- вертикальная фильтрация ----------
#pragma omp parallel for default(none) shared(n, m, temp, out)
  for (int yi = 0; yi < m; ++yi) {
    for (int xy = 0; xy < n; ++xy) {
      out[(yi * n) + xy] = VerticalFilterAtOMP(temp, n, m, xy, yi);
    }
  }

#pragma omp barrier

#if _OPENMP >= 201811
  omp_pause_resource_all(omp_pause_hard);
#endif

  return true;
}

bool ChaschinVLinearFiltrationOMP::PostProcessingImpl() {
  return true;
}

}  // namespace chaschin_v_linear_image_filtration_omp
