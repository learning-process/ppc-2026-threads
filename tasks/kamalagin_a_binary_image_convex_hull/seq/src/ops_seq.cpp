#include "kamalagin_a_binary_image_convex_hull/seq/include/ops_seq.hpp"

#include <cstddef>

#include "kamalagin_a_binary_image_convex_hull/common/include/common.hpp"
#include "kamalagin_a_binary_image_convex_hull/common/include/run_impl.hpp"

namespace kamalagin_a_binary_image_convex_hull {

KamalaginABinaryImageConvexHullSEQ::KamalaginABinaryImageConvexHullSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = HullList{};
}

bool KamalaginABinaryImageConvexHullSEQ::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows < 0 || img.cols < 0) {
    return false;
  }
  if (img.rows == 0 || img.cols == 0) {
    return img.data.empty();
  }
  if (img.rows > 1000 || img.cols > 1000) {
    return false;
  }
  return (static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols)) == img.data.size();
}

bool KamalaginABinaryImageConvexHullSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool KamalaginABinaryImageConvexHullSEQ::RunImpl() {
  detail::RunBinaryImageConvexHull(GetInput(), GetOutput());
  return true;
}

bool KamalaginABinaryImageConvexHullSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kamalagin_a_binary_image_convex_hull
