#include "dergynov_s_integrals_multistep_rectangle/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace dergynov_s_integrals_multistep_rectangle {

DergynovSIntegralsMultistepRectangleSEQ::DergynovSIntegralsMultistepRectangleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleSEQ::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();
  
  if (!func) return false;
  
  if (n <= 0) return false;
  
  if (borders.empty()) return false;
  
  for (const auto &[left, right] : borders) {
    if (!std::isfinite(left) || !std::isfinite(right)) return false;
    if (left >= right) return false;
  }
  
  return true;
}

bool DergynovSIntegralsMultistepRectangleSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

namespace {
bool NextIndex(std::vector<int> &idx, int dim, int n) {
  for (int pos = 0; pos < dim; ++pos) {
    ++idx[pos];
    if (idx[pos] < n) return true;
    idx[pos] = 0;
  }
  return false;
}
}  // namespace

bool DergynovSIntegralsMultistepRectangleSEQ::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  
  for (int i = 0; i < dim; ++i) {
    const double left = borders[i].first;
    const double right = borders[i].second;
    h[i] = (right - left) / n;
    cell_volume *= h[i];
  }

  std::vector<int> idx(dim, 0);
  std::vector<double> point(dim);
  double sum = 0.0;

  do {
    for (int i = 0; i < dim; ++i) {
      point[i] = borders[i].first + (idx[i] + 0.5) * h[i];
    }
    
    double f_val = func(point);
    if (!std::isfinite(f_val)) return false;
    
    sum += f_val;
    
  } while (NextIndex(idx, dim, n));

  GetOutput() = sum * cell_volume;
  
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleSEQ::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle