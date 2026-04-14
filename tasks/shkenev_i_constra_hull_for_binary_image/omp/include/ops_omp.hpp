#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "shkenev_i_constra_hull_for_binary_image/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkenev_i_constra_hull_for_binary_image {

class ShkenevIConstrHullOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ShkenevIConstrHullOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ThresholdImage();
  void FindComponents();
  static std::vector<Point> BuildHull(const std::vector<Point> &points);
  static size_t Index(int x, int y, int width);
  void ExploreComponent(int start_col, int start_row, int width, int height, std::vector<uint8_t> &visited,
                        std::vector<Point> &component);

  BinaryImage work_;
};

}  // namespace shkenev_i_constra_hull_for_binary_image
