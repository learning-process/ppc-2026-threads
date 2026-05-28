#include "shakirova_e_sobel_edge_detection/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "shakirova_e_sobel_edge_detection/common/include/common.hpp"

namespace shakirova_e_sobel_edge_detection {

ShakirovaESobelEdgeDetectionTBB::ShakirovaESobelEdgeDetectionTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ShakirovaESobelEdgeDetectionTBB::ValidationImpl() {
  return GetInput().IsValid();
}

bool ShakirovaESobelEdgeDetectionTBB::PreProcessingImpl() {
  const auto &img = GetInput();
  GetOutput().assign(static_cast<size_t>(img.width) * static_cast<size_t>(img.height), 0);
  return true;
}

bool ShakirovaESobelEdgeDetectionTBB::RunImpl() {
  const auto &img = GetInput();
  const int h = img.height;
  const int w = img.width;
  const int *inp = img.pixels.data();
  int *out = GetOutput().data();

  tbb::parallel_for(tbb::blocked_range<int>(1, h - 1, 8), [inp, out, w](const tbb::blocked_range<int> &r) {
    for (int row = r.begin(); row < r.end(); ++row) {
      const int *prev = inp + static_cast<ptrdiff_t>((row - 1) * w);
      const int *curr = inp + static_cast<ptrdiff_t>((row)*w);
      const int *next = inp + static_cast<ptrdiff_t>((row + 1) * w);
      int *out_row = out + static_cast<ptrdiff_t>((row)*w);

      for (int col = 1; col < w - 1; ++col) {
        const int gx =
            -prev[col - 1] + prev[col + 1] - (2 * curr[col - 1]) + (2 * curr[col + 1]) - next[col - 1] + next[col + 1];
        const int gy =
            -prev[col - 1] - (2 * prev[col]) - prev[col + 1] + next[col - 1] + (2 * next[col]) + next[col + 1];

        const int agx = std::abs(gx);
        const int agy = std::abs(gy);
        const int magnitude = ((std::max(agx, agy) * 123) + (std::min(agx, agy) * 51)) >> 7;
        out_row[col] = magnitude > 255 ? 255 : magnitude;
      }
    }
  });

  return true;
}

bool ShakirovaESobelEdgeDetectionTBB::PostProcessingImpl() {
  return true;
}

}  // namespace shakirova_e_sobel_edge_detection
