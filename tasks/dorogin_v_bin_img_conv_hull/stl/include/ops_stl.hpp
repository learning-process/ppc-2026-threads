#pragma once

#include <cstddef>
#include <vector>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull {

class DoroginVBinImgConvHullSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit DoroginVBinImgConvHullSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ThresholdImage();
  void FindComponents();
  static void FillConvexHullsChunk(const std::vector<std::vector<Point>> &components,
                                   std::vector<std::vector<Point>> &convex_hulls, std::size_t begin, std::size_t end);
  static std::vector<Point> BuildHull(const std::vector<Point> &points);
  static std::size_t Index(int col, int row, int width);
  void ExploreComponent(int start_col, int start_row, int width, int height, std::vector<bool> &visited,
                        std::vector<Point> &component);

  BinaryImage w_;
};

}  // namespace dorogin_v_bin_img_conv_hull
