#include "korolev_k_matrix_mult/stl/include/ops_stl.hpp"

#include <cstddef>
#include <functional>
#include <thread>
#include <vector>

#include "korolev_k_matrix_mult/common/include/common.hpp"
#include "korolev_k_matrix_mult/common/include/strassen_impl.hpp"

namespace korolev_k_matrix_mult {

namespace {

void ParallelInvokeThreads(std::vector<std::function<void()>> &tasks) {
  std::vector<std::thread> threads;
  threads.reserve(tasks.size());
  for (auto &t : tasks) {
    threads.emplace_back(t);
  }
  for (auto &th : threads) {
    th.join();
  }
}

void PadInputBlocks(const InType &in, size_t n, size_t np2, std::vector<double> &a_pad, std::vector<double> &b_pad) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      a_pad[(i * np2) + j] = in.A[(i * n) + j];
      b_pad[(i * np2) + j] = in.B[(i * n) + j];
    }
  }
}

void CopyResultCorner(const std::vector<double> &c_pad, std::vector<double> &out, size_t n, size_t np2) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      out[(i * n) + j] = c_pad[(i * np2) + j];
    }
  }
}

}  // namespace

KorolevKMatrixMultSTL::KorolevKMatrixMultSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KorolevKMatrixMultSTL::ValidationImpl() {
  const auto &in = GetInput();
  return in.n > 0 && in.A.size() == in.n * in.n && in.B.size() == in.n * in.n && GetOutput().empty();
}

bool KorolevKMatrixMultSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().n * GetInput().n);
  return true;
}

bool KorolevKMatrixMultSTL::RunImpl() {
  const auto &in = GetInput();
  size_t n = in.n;
  size_t np2 = strassen_impl::NextPowerOf2(n);

  if (np2 == n) {
    strassen_impl::StrassenMultiply(in.A, in.B, GetOutput(), n, ParallelInvokeThreads);
  } else {
    std::vector<double> a_pad(np2 * np2, 0);
    std::vector<double> b_pad(np2 * np2, 0);
    std::vector<double> c_pad(np2 * np2, 0);
    PadInputBlocks(in, n, np2, a_pad, b_pad);
    strassen_impl::StrassenMultiply(a_pad, b_pad, c_pad, np2, ParallelInvokeThreads);
    CopyResultCorner(c_pad, GetOutput(), n, np2);
  }
  return true;
}

bool KorolevKMatrixMultSTL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace korolev_k_matrix_mult
