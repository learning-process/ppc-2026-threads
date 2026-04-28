#include "korolev_k_matrix_mult/omp/include/ops_omp.hpp"

#include <cstddef>
#include <functional>
#include <vector>

#include "korolev_k_matrix_mult/common/include/common.hpp"
#include "korolev_k_matrix_mult/common/include/strassen_impl.hpp"
#include "util/include/util.hpp"

namespace korolev_k_matrix_mult {

KorolevKMatrixMultOMP::KorolevKMatrixMultOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KorolevKMatrixMultOMP::ValidationImpl() {
  const auto &in = GetInput();
  return in.n > 0 && in.A.size() == in.n * in.n && in.B.size() == in.n * in.n && GetOutput().empty();
}

bool KorolevKMatrixMultOMP::PreProcessingImpl() {
  GetOutput().resize(GetInput().n * GetInput().n);
  return true;
}

bool KorolevKMatrixMultOMP::RunImpl() {
  const auto &in = GetInput();
  size_t n = in.n;
  size_t np2 = strassen_impl::NextPowerOf2(n);

  auto parallel_run = [](std::vector<std::function<void()>> &tasks) {
#pragma omp parallel default(none) num_threads(ppc::util::GetNumThreads()) shared(tasks)
#pragma omp single
    {
      for (size_t i = 0; i < tasks.size(); ++i) {
#pragma omp task default(none) firstprivate(i) shared(tasks)
        tasks[i]();
      }
#pragma omp taskwait
    }
  };

  if (np2 == n) {
    strassen_impl::StrassenMultiply(in.A, in.B, GetOutput(), n, parallel_run);
  } else {
    std::vector<double> a_pad(np2 * np2, 0);
    std::vector<double> b_pad(np2 * np2, 0);
    std::vector<double> c_pad(np2 * np2, 0);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        a_pad[(i * np2) + j] = in.A[(i * n) + j];
        b_pad[(i * np2) + j] = in.B[(i * n) + j];
      }
    }
    strassen_impl::StrassenMultiply(a_pad, b_pad, c_pad, np2, parallel_run);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        GetOutput()[(i * n) + j] = c_pad[(i * np2) + j];
      }
    }
  }
  return true;
}

bool KorolevKMatrixMultOMP::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace korolev_k_matrix_mult
