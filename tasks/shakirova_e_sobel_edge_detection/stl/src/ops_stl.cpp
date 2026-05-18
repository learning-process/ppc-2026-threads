#include "shakirova_e_sobel_edge_detection/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "shakirova_e_sobel_edge_detection/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shakirova_e_sobel_edge_detection {

ShakirovaESobelEdgeDetectionSTL::ShakirovaESobelEdgeDetectionSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ShakirovaESobelEdgeDetectionSTL::ValidationImpl() { return GetInput().IsValid(); }

bool ShakirovaESobelEdgeDetectionSTL::PreProcessingImpl() {
  const auto &img = GetInput();
  GetOutput().assign(static_cast<size_t>(img.width) * static_cast<size_t>(img.height), 0);
  return true;
}

bool ShakirovaESobelEdgeDetectionSTL::RunImpl() {
  const auto &img = GetInput();
  const int h = img.height;
  const int w = img.width;
  const int *inp = img.pixels.data();
  int *out = GetOutput().data();

  const int inner_rows = h - 2;
  if (inner_rows <= 0) {
    return true;
  }

  const int num_t = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_t));

  for (int tid = 0; tid < num_t; ++tid) {
    threads.emplace_back([inp, out, w, inner_rows, num_t, tid]() {
      const int chunk = inner_rows / num_t;
      const int row_begin = tid * chunk + 1;
      const int row_end = (tid == num_t - 1) ? (inner_rows + 1) : (row_begin + chunk);

      for (int row = row_begin; row < row_end; ++row) {
        const int *prev = inp + static_cast<ptrdiff_t>((row - 1) * w);
        const int *curr = inp + static_cast<ptrdiff_t>((row) * w);
        const int *next = inp + static_cast<ptrdiff_t>((row + 1) * w);
        int *out_row    = out + static_cast<ptrdiff_t>((row) * w);

        for (int col = 1; col < w - 1; ++col) {
          const int gx = -prev[col - 1] + prev[col + 1]
                         - (2 * curr[col - 1]) + (2 * curr[col + 1])
                         - next[col - 1] + next[col + 1];
          const int gy = -prev[col - 1] - (2 * prev[col]) - prev[col + 1]
                         + next[col - 1] + (2 * next[col]) + next[col + 1];

          const int agx = std::abs(gx);
          const int agy = std::abs(gy);
          const int magnitude = ((std::max(agx, agy) * 123) + (std::min(agx, agy) * 51)) >> 7;
          out_row[col] = magnitude > 255 ? 255 : magnitude;
        }
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  return true;
}

bool ShakirovaESobelEdgeDetectionSTL::PostProcessingImpl() { return true; }

}  // namespace shakirova_e_sobel_edge_detection 
