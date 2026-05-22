#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull {

class DoroginVBinImgConvHullALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit DoroginVBinImgConvHullALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ThresholdPixelsParallel();
  void CollectComponentsSequential();
  static std::vector<Point> BuildHull(const std::vector<Point> &points);
  static std::size_t Index(int col, int row, int width);
  void FloodFill(int seed_x, int seed_y, int width, int height, std::vector<std::uint8_t> &visited,
                 std::vector<Point> &component);

  BinaryImage work_;
};

}  // namespace dorogin_v_bin_img_conv_hull
