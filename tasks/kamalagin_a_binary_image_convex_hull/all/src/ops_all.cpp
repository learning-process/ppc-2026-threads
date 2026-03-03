#include "kamalagin_a_binary_image_convex_hull/all/include/ops_all.hpp"

#include "kamalagin_a_binary_image_convex_hull/common/include/run_impl.hpp"

namespace kamalagin_a_binary_image_convex_hull {

KamalaginABinaryImageConvexHullALL::KamalaginABinaryImageConvexHullALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = HullList{};
}

bool KamalaginABinaryImageConvexHullALL::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows <= 0 || img.cols <= 0) {
    return false;
  }
  if (img.rows > 1000 || img.cols > 1000) {
    return false;
  }
  return static_cast<size_t>(img.rows * img.cols) == img.data.size();
}

bool KamalaginABinaryImageConvexHullALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool KamalaginABinaryImageConvexHullALL::RunImpl() {
  detail::RunBinaryImageConvexHull(GetInput(), GetOutput());
  return true;
}

bool KamalaginABinaryImageConvexHullALL::PostProcessingImpl() {
  return true;
}

}  // namespace kamalagin_a_binary_image_convex_hull
