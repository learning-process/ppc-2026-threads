#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "paramonov_v_bin_img_conv_hul_omp/common/include/common.hpp"
#include "task/include/task.hpp"

namespace paramonov_v_bin_img_conv_hul_omp {

class ConvexHullOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit ConvexHullOMP(const InputType &input);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  void BinarizeImage(uint8_t threshold = 128);
  void ExtractConnectedComponents();

  [[nodiscard]] static std::vector<PixelPoint> ComputeConvexHull(const std::vector<PixelPoint> &points);
  [[nodiscard]] static int64_t Orientation(const PixelPoint &p, const PixelPoint &q, const PixelPoint &r);
  static size_t PixelIndex(int row, int col, int cols) {
    return (static_cast<size_t>(row) * static_cast<size_t>(cols)) + static_cast<size_t>(col);
  }

  void FloodFill(int start_row, int start_col, std::vector<bool> &visited, std::vector<PixelPoint> &component) const;

  InputType working_image_;
};

}  // namespace paramonov_v_bin_img_conv_hul_omp
