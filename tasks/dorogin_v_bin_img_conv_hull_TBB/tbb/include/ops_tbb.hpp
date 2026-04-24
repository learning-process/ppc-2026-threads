#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_TBB/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull_tbb {

class DoroginVImgConvHullTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit DoroginVImgConvHullTBB(const InputType &input);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static size_t PixelIndex(int row, int col, int cols) {
    return (static_cast<size_t>(row) * static_cast<size_t>(cols)) + static_cast<size_t>(col);
  }

  static int64_t Cross(const PixelPoint &a, const PixelPoint &b, const PixelPoint &c);
  static std::vector<PixelPoint> ComputeConvexHull(const std::vector<PixelPoint> &points);

  void BinarizeImage(uint8_t threshold = 127);
  void FloodFillComponent(int row_start, int col_start, std::vector<bool> *visited,
                          std::vector<PixelPoint> *component) const;
  void ExtractComponentsAndHulls();

  InputType image_;
};

}  // namespace dorogin_v_bin_img_conv_hull_tbb
