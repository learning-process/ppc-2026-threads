#include "kamalagin_a_binary_image_convex_hull_omp/omp/include/ops_omp.hpp"

#include <cstddef>

#include "kamalagin_a_binary_image_convex_hull_omp/common/include/common.hpp"
#include "kamalagin_a_binary_image_convex_hull_omp/common/include/run_impl.hpp"

namespace kamalagin_a_binary_image_convex_hull_omp {

KamalaginABinaryImageConvexHullOMP::KamalaginABinaryImageConvexHullOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KamalaginABinaryImageConvexHullOMP::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows < 0 || img.cols < 0) {
    return false;
  }
  if (img.rows == 0 || img.cols == 0) {
    return img.data.empty();
  }
  return (static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols)) == img.data.size();
}

bool KamalaginABinaryImageConvexHullOMP::PreProcessingImpl() {
  return true;
}

bool KamalaginABinaryImageConvexHullOMP::RunImpl() {
  RunBinaryImageConvexHullOmp(GetInput(), GetOutput());
  return true;
}

bool KamalaginABinaryImageConvexHullOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kamalagin_a_binary_image_convex_hull_omp
