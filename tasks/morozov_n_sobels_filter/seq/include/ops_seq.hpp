#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

#include "morozov_n_sobels_filter/common/include/common.hpp"
#include "task/include/task.hpp"

namespace morozov_n_sobels_filter {

class MorozovNSobelsFilterSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit MorozovNSobelsFilterSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void Filter(const Image &img);
  static uint8_t CalculateNewPixelColor(const Image &img, size_t x, size_t y);

  static constexpr std::array<std::array<int, 3>, 3> kKernelX = {
      std::array<int, 3>{-1, 0, 1}, std::array<int, 3>{-2, 0, 2}, std::array<int, 3>{-1, 0, 1}};

  static constexpr std::array<std::array<int, 3>, 3> kKernelY = {
      std::array<int, 3>{-1, -2, -1}, std::array<int, 3>{0, 0, 0}, std::array<int, 3>{1, 2, 1}};
  Image result_image_;
};

}  // namespace morozov_n_sobels_filter
