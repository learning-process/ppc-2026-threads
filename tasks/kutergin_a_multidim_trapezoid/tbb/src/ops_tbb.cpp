#include "kutergin_a_multidim_trapezoid/tbb/include/ops_tbb.hpp"

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace kutergin_a_multidim_trapezoid {
namespace {

bool ValidateBorders(const std::vector<std::pair<double, double>> &borders) {
  return std::ranges::all_of(
      borders, [](const auto &p) { return std::isfinite(p.first) && std::isfinite(p.second) && (p.first < p.second); });
}

}  // namespace

KuterginAMultidimTrapezoidTBB::KuterginAMultidimTrapezoidTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool KuterginAMultidimTrapezoidTBB::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();
  return func && n > 0 && !borders.empty() && ValidateBorders(borders);
}

bool KuterginAMultidimTrapezoidTBB::PreProcessingImpl() {
    local_input_ = GetInput(); // Копируем всё: функцию, границы, n
    res_ = 0.0;
    return true;
}

bool KuterginAMultidimTrapezoidTBB::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  for (int i = 0; i < dim; ++i) {
    h[i] = (borders[i].second - borders[i].first) / n;
    cell_volume *= h[i];
  }

  size_t total_points = 1;
  for (int i = 0; i < dim; ++i) total_points *= (static_cast<size_t>(n) + 1);

  double total_sum = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, total_points), 
      0.0,
      [=](const tbb::blocked_range<size_t> &r, double local_sum) {
        std::vector<double> point(dim); 
        
        for (size_t i = r.begin(); i < r.end(); ++i) {
          double weight = 1.0;
          size_t temp_idx = i;

          for (int d = 0; d < dim; ++d) {
            int coord_idx = static_cast<int>(temp_idx % (n + 1));
            temp_idx /= (n + 1);

            point[d] = borders[d].first + coord_idx * h[d];
            
            // Расчет веса для метода трапеций по каждому измерению
            if (coord_idx == 0 || coord_idx == n) {
              weight *= 0.5;
            }
          }

          local_sum += weight * func(point);
        }
        return local_sum;
      },
      std::plus<double>()
  );

  GetOutput() = total_sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool KuterginAMultidimTrapezoidTBB::PostProcessingImpl() {
  GetOutput() = res_; // Только здесь отдаем результат в систему
  return true;
}

}  // namespace kutergin_a_multidim_trapezoid