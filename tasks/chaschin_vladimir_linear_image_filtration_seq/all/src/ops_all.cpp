#include "chaschin_vladimir_linear_image_filtration_seq/all/include/ops_all.hpp"

#include <omp.h>

#include <thread>
#include <utility>
#include <vector>

namespace chaschin_v_linear_image_filtration_all {

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
  out.resize(static_cast<std::vector<float>::size_type>(n) * static_cast<std::vector<float>::size_type>(m));

  std::vector<float> temp(static_cast<std::vector<float>::size_type>(n) *
                          static_cast<std::vector<float>::size_type>(m));

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 4;
  }
  if (num_threads > static_cast<unsigned int>(m)) {
    num_threads = m;
  }

  auto worker_horizontal = [&](int start_y, int end_y) {
    for (int yi = start_y; yi < end_y; ++yi) {
#pragma omp simd
      for (int xf = 0; xf < n; ++xf) {
        temp[(yi * n) + xf] = HorizontalFilterAtSTL(image, n, xf, yi);
      }
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  int chunk = m / num_threads;
  int remainder = m % num_threads;
  int current_y = 0;

  for (unsigned int i = 0; i < num_threads; ++i) {
    int end_y = current_y + chunk + (static_cast<int>(i) < remainder ? 1 : 0);
    threads.emplace_back(worker_horizontal, current_y, end_y);
    current_y = end_y;
  }
  for (auto &t : threads) {
    t.join();
  }

  auto worker_vertical = [&](int start_y, int end_y) {
    for (int yi = start_y; yi < end_y; ++yi) {
#pragma omp simd
      for (int xy = 0; xy < n; ++xy) {
        out[(yi * n) + xy] = VerticalFilterAtSTL(temp, n, m, xy, yi);
      }
    }
  };

  threads.clear();
  current_y = 0;

  for (unsigned int i = 0; i < num_threads; ++i) {
    int end_y = current_y + chunk + (static_cast<int>(i) < remainder ? 1 : 0);
    threads.emplace_back(worker_vertical, current_y, end_y);
    current_y = end_y;
  }
  for (auto &t : threads) {
    t.join();
  }

  return true;
}

bool ChaschinVLinearFiltrationSTL::PostProcessingImpl() {
  return true;
}

}  // namespace chaschin_v_linear_image_filtration_all
