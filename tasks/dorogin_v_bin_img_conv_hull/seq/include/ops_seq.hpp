#pragma once

#include <cstddef>
#include <vector>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull {

class DoroginVBinImgConvHullSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit DoroginVBinImgConvHullSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ThresholdImage();
  void FindComponents();
  static std::vector<Point> BuildHull(const std::vector<Point> &points);
  static size_t Index(int x, int y, int width);
  void ExploreComponent(int start_col, int start_row, int width, int height, std::vector<bool> &visited,
                        std::vector<Point> &component);

  BinaryImage w_;
};

}  // namespace dorogin_v_bin_img_conv_hull
