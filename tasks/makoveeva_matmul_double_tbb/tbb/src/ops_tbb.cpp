#include "makoveeva_matmul_double_tbb/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_tbb/common/include/common.hpp"

namespace makoveeva_matmul_double_tbb {

MatmulDoubleTBBTask::MatmulDoubleTBBTask(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool MatmulDoubleTBBTask::ValidationImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);

  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool MatmulDoubleTBBTask::PreProcessingImpl() {
  const auto &input = GetInput();
  n_ = std::get<0>(input);
  A_ = std::get<1>(input);
  B_ = std::get<2>(input);
  C_.assign(n_ * n_, 0.0);

  return true;
}

bool MatmulDoubleTBBTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const size_t n = n_;
  const auto &a = A_;
  const auto &b = B_;
  auto &c = C_;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      for (size_t j = 0; j < n; ++j) {
        double sum = 0.0;
        for (size_t k = 0; k < n; ++k) {
          sum += a[(i * n) + k] * b[(k * n) + j];
        }
        c[(i * n) + j] = sum;
      }
    }
  });

  this->GetOutput() = C_;
  return true;
}

bool MatmulDoubleTBBTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_tbb
