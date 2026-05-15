#include "chaschin_v_linear_image_filtration_stl/stl/include/ops_stl.hpp"

#include <algorithm>
#include <future>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

#include "chaschin_v_linear_image_filtration_seq/common/include/common.hpp"

namespace chaschin_v_linear_image_filtration_stl {

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

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) num_threads = 2;

  // ---------- горизонтальная фильтрация ----------
  {
    std::vector<std::future<void>> futures;
    for (unsigned int t = 0; t < num_threads; ++t) {
      futures.push_back(std::async(std::launch::async, [t, num_threads, m, n, &image, &temp]() {
        for (int yi = t; yi < m; yi += num_threads) {
          for (int xf = 0; xf < n; ++xf) {
            temp[(yi * n) + xf] = HorizontalFilterAtSTL(image, n, xf, yi);
          }
        }
      }));
    }
    for (auto &f : futures) {
      f.wait();
    }
  }

  // ---------- вертикальная фильтрация ----------
  {
    std::vector<std::future<void>> futures;
    for (unsigned int t = 0; t < num_threads; ++t) {
      futures.push_back(std::async(std::launch::async, [t, num_threads, m, n, &temp, &out]() {
        for (int yi = t; yi < m; yi += num_threads) {
          for (int xy = 0; xy < n; ++xy) {
            out[(yi * n) + xy] = VerticalFilterAtSTL(temp, n, m, xy, yi);
          }
        }
      }));
    }
    for (auto &f : futures) {
      f.wait();
    }
  }

  return true;
}

bool ChaschinVLinearFiltrationSTL::PostProcessingImpl() {
  return true;
}

}  // namespace chaschin_v_linear_image_filtration_stl