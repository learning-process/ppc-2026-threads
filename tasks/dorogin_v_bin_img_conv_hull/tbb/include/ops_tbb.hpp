#pragma once

#include <cstddef>
#include <vector>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull {

class DoroginVBinImgConvHullTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit DoroginVBinImgConvHullTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ThresholdImage();
  void FindComponents();
  static std::vector<Point> BuildHull(const std::vector<Point> &points);
  static std::size_t Index(int col, int row, int width);
  void ExploreComponent(int start_col, int start_row, int width, int height, std::vector<bool> &visited,
                        std::vector<Point> &component);

  BinaryImage w_;
};

}  // namespace dorogin_v_bin_img_conv_hull
