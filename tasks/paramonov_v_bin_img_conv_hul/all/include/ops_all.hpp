#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "paramonov_v_bin_img_conv_hul/all/include/common.hpp"
#include "task/include/task.hpp"

namespace paramonov_v_bin_img_conv_hul_all {

class ConvexHullAll : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    // Пробуем разные варианты:
    // return ppc::task::TypeOfTask::kAll;      // если есть
    // return ppc::task::TypeOfTask::kGeneric;  // если есть
    return ppc::task::TypeOfTask::kSTL;  // временно, чтобы собралось
  }

  explicit ConvexHullAll(const InputType &input);

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

  [[nodiscard]] static PixelPoint FindLowestPoint(const std::vector<PixelPoint> &points);
  [[nodiscard]] static std::vector<PixelPoint> SortPointsByAngle(const std::vector<PixelPoint> &points,
                                                                 const PixelPoint &lowest_point);
  [[nodiscard]] static std::vector<PixelPoint> RemoveCollinearPoints(const std::vector<PixelPoint> &sorted_points,
                                                                     const PixelPoint &lowest_point);
  [[nodiscard]] static std::vector<PixelPoint> BuildHull(const std::vector<PixelPoint> &unique_points,
                                                         const PixelPoint &lowest_point);
  [[nodiscard]] static std::vector<PixelPoint> HandleCollinearCase(const std::vector<PixelPoint> &points,
                                                                   const PixelPoint &lowest_point);

  InputType working_image_;
};

}  // namespace paramonov_v_bin_img_conv_hul_all
