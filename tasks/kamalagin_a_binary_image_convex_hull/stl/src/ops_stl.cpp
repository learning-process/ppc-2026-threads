#include "kamalagin_a_binary_image_convex_hull/stl/include/ops_stl.hpp"

#include "kamalagin_a_binary_image_convex_hull/common/include/run_impl.hpp"

namespace kamalagin_a_binary_image_convex_hull {

KamalaginABinaryImageConvexHullSTL::KamalaginABinaryImageConvexHullSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = HullList{};
}

bool KamalaginABinaryImageConvexHullSTL::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows <= 0 || img.cols <= 0) {
    return false;
  }
  if (img.rows > 1000 || img.cols > 1000) {
    return false;
  }
  return static_cast<size_t>(img.rows * img.cols) == img.data.size();
}

bool KamalaginABinaryImageConvexHullSTL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool KamalaginABinaryImageConvexHullSTL::RunImpl() {
  detail::RunBinaryImageConvexHull(GetInput(), GetOutput());
  return true;
}

bool KamalaginABinaryImageConvexHullSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kamalagin_a_binary_image_convex_hull
