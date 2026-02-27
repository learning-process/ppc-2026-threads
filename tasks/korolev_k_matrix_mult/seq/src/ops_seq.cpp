#include "korolev_k_matrix_mult/seq/include/ops_seq.hpp"

#include <cstddef>
#include <vector>

#include "korolev_k_matrix_mult/common/include/common.hpp"
#include "korolev_k_matrix_mult/common/include/strassen_impl.hpp"

namespace korolev_k_matrix_mult {

KorolevKMatrixMultSEQ::KorolevKMatrixMultSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KorolevKMatrixMultSEQ::ValidationImpl() {
  const auto &in = GetInput();
  return in.n > 0 && in.A.size() == in.n * in.n && in.B.size() == in.n * in.n && GetOutput().empty();
}

bool KorolevKMatrixMultSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().n * GetInput().n);
  return true;
}

bool KorolevKMatrixMultSEQ::RunImpl() {
  const auto &in = GetInput();
  size_t n = in.n;
  size_t np2 = strassen_impl::NextPowerOf2(n);

  auto sequential_run = [](std::vector<std::function<void()>> &tasks) {
    for (auto &t : tasks) t();
  };

  if (np2 == n) {
    strassen_impl::StrassenMultiply(in.A, in.B, GetOutput(), n, sequential_run);
  } else {
    std::vector<double> A_pad(np2 * np2, 0), B_pad(np2 * np2, 0), C_pad(np2 * np2, 0);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        A_pad[i * np2 + j] = in.A[i * n + j];
        B_pad[i * np2 + j] = in.B[i * n + j];
      }
    }
    strassen_impl::StrassenMultiply(A_pad, B_pad, C_pad, np2, sequential_run);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        GetOutput()[i * n + j] = C_pad[i * np2 + j];
      }
    }
  }
  return true;
}

bool KorolevKMatrixMultSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace korolev_k_matrix_mult
