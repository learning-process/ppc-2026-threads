#include "dorogin_v_bin_img_conv_hull_TBB/tbb/include/ops_tbb.hpp"

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <stack>
#include <utility>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_TBB/common/include/common.hpp"

namespace dorogin_v_bin_img_conv_hull_tbb {
namespace {
constexpr std::array<std::pair<int, int>, 4> k4Connected = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};
}  // namespace

DoroginVImgConvHullTBB::DoroginVImgConvHullTBB(const InputType &input) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input;
}

bool DoroginVImgConvHullTBB::ValidationImpl() {
  const auto &input = GetInput();
  if (input.rows <= 0 || input.cols <= 0) {
    return false;
  }
  const size_t required_size = static_cast<size_t>(input.rows) * static_cast<size_t>(input.cols);
  return input.pixels.size() == required_size;
}

bool DoroginVImgConvHullTBB::PreProcessingImpl() {
  image_ = GetInput();
  BinarizeImage();
  GetOutput().clear();
  return true;
}

bool DoroginVImgConvHullTBB::RunImpl() {
  ExtractComponentsAndHulls();
  return true;
}

bool DoroginVImgConvHullTBB::PostProcessingImpl() {
  return true;
}

void DoroginVImgConvHullTBB::BinarizeImage(uint8_t threshold) {
  auto &pixels = image_.pixels;
  const size_t pixel_count = pixels.size();
  tbb::parallel_for(size_t{0}, pixel_count,
                    [&](size_t i) { pixels[i] = pixels[i] > threshold ? uint8_t{255} : uint8_t{0}; });
}

void DoroginVImgConvHullTBB::FloodFillComponent(int row_start, int col_start, std::vector<bool> *visited,
                                                std::vector<PixelPoint> *component) const {
  std::stack<PixelPoint> stack;
  stack.emplace(row_start, col_start);
  (*visited)[PixelIndex(row_start, col_start, image_.cols)] = true;

  while (!stack.empty()) {
    PixelPoint cur = stack.top();
    stack.pop();
    component->push_back(cur);

    for (const auto &[dr, dc] : k4Connected) {
      const int nr = cur.row + dr;
      const int nc = cur.col + dc;
      if (nr < 0 || nr >= image_.rows || nc < 0 || nc >= image_.cols) {
        continue;
      }

      const size_t idx = PixelIndex(nr, nc, image_.cols);
      if (!(*visited)[idx] && image_.pixels[idx] == 255) {
        (*visited)[idx] = true;
        stack.emplace(nr, nc);
      }
    }
  }
}

int64_t DoroginVImgConvHullTBB::Cross(const PixelPoint &a, const PixelPoint &b, const PixelPoint &c) {
  const int64_t abx = static_cast<int64_t>(b.col) - static_cast<int64_t>(a.col);
  const int64_t aby = static_cast<int64_t>(b.row) - static_cast<int64_t>(a.row);
  const int64_t acx = static_cast<int64_t>(c.col) - static_cast<int64_t>(a.col);
  const int64_t acy = static_cast<int64_t>(c.row) - static_cast<int64_t>(a.row);
  return (abx * acy) - (aby * acx);
}

std::vector<PixelPoint> DoroginVImgConvHullTBB::ComputeConvexHull(const std::vector<PixelPoint> &points) {
  if (points.size() <= 2) {
    return points;
  }

  std::vector<PixelPoint> ordered = points;
  tbb::parallel_sort(ordered.begin(), ordered.end(), [](const PixelPoint &a, const PixelPoint &b) {
    if (a.row != b.row) {
      return a.row < b.row;
    }
    return a.col < b.col;
  });
  const auto unique_end = std::ranges::unique(
      ordered, [](const PixelPoint &a, const PixelPoint &b) { return a.row == b.row && a.col == b.col; });
  ordered.erase(unique_end.begin(), unique_end.end());
  if (ordered.size() <= 2) {
    return ordered;
  }

  std::vector<PixelPoint> lower;
  for (const auto &p : ordered) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower.back(), p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  std::vector<PixelPoint> upper;
  for (const auto &point : std::ranges::reverse_view(ordered)) {
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper.back(), point) <= 0) {
      upper.pop_back();
    }
    upper.push_back(point);
  }

  lower.pop_back();
  upper.pop_back();
  lower.insert(lower.end(), upper.begin(), upper.end());
  return lower;
}

void DoroginVImgConvHullTBB::ExtractComponentsAndHulls() {
  const size_t total = static_cast<size_t>(image_.rows) * static_cast<size_t>(image_.cols);
  std::vector<bool> visited(total, false);
  auto &hulls = GetOutput();

  for (int row_idx = 0; row_idx < image_.rows; ++row_idx) {
    for (int col_idx = 0; col_idx < image_.cols; ++col_idx) {
      const size_t idx = PixelIndex(row_idx, col_idx, image_.cols);
      if (image_.pixels[idx] != 255 || visited[idx]) {
        continue;
      }

      std::vector<PixelPoint> component;
      FloodFillComponent(row_idx, col_idx, &visited, &component);
      if (!component.empty()) {
        hulls.push_back(ComputeConvexHull(component));
      }
    }
  }
}

}  // namespace dorogin_v_bin_img_conv_hull_tbb
