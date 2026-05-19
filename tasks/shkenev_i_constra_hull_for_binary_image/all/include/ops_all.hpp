#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "shkenev_i_constra_hull_for_binary_image/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkenev_i_constra_hull_for_binary_image {

class ShkenevIConstrHullALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit ShkenevIConstrHullALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ThresholdImage();
  void FindComponentsTBB();

  static std::vector<Point> BuildHull(const std::vector<Point> &points);
  static size_t Index(int x, int y, int width);

  void BFS(int sx, int sy, int width, int height, std::vector<uint8_t> &visited, std::vector<Point> &comp);

  BinaryImage work_;
};

}  // namespace shkenev_i_constra_hull_for_binary_image
