#pragma once

#include <cstddef>  // для size_t
#include <cstdint>  // для uint8_t, int64_t
#include <vector>   // для std::vector

#include "paramonov_v_bin_img_conv_hul/common/include/common.hpp"
#include "task/include/task.hpp"

namespace paramonov_v_bin_img_conv_hul {

class ConvexHullSequential : public HullTaskBase {
 public:
  static constexpr ppc::task::TypeOfTask StaticTaskType() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit ConvexHullSequential(const InputType &input);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void BinarizeImage(uint8_t threshold = 128);
  void ExtractConnectedComponents();
  [[nodiscard]] std::vector<PixelPoint> ComputeConvexHull(const std::vector<PixelPoint> &points) const;
  void FloodFill(int start_row, int start_col, std::vector<bool> &visited, std::vector<PixelPoint> &component) const;

  [[nodiscard]] static int64_t Orientation(const PixelPoint &p, const PixelPoint &q, const PixelPoint &r);
  static size_t PixelIndex(int row, int col, int cols) {
    return (static_cast<size_t>(row) * static_cast<size_t>(cols)) + static_cast<size_t>(col);
  }

  InputType working_image_;
};

}  // namespace paramonov_v_bin_img_conv_hul
