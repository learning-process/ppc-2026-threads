#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "morozov_n_sobels_filter/common/include/common.hpp"
#include "task/include/task.hpp"

namespace morozov_n_sobels_filter {

class MorozovNSobelsFilterALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit MorozovNSobelsFilterALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void SplitRows(size_t global_rows, size_t proc_num, size_t &start, size_t &count) const;
  void SendImageDataFromZeroProc(const Image &global, size_t halo);
  void CollectResult();

  static void Filter(const Image &img, Image &local_result, size_t start_row, size_t end_row);
  static uint8_t CalculateNewPixelColor(const Image &img, size_t x, size_t y);

  static constexpr std::array<std::array<int, 3>, 3> kKernelX = {
      std::array<int, 3>{-1, 0, 1}, std::array<int, 3>{-2, 0, 2}, std::array<int, 3>{-1, 0, 1}};

  static constexpr std::array<std::array<int, 3>, 3> kKernelY = {
      std::array<int, 3>{-1, -2, -1}, std::array<int, 3>{0, 0, 0}, std::array<int, 3>{1, 2, 1}};

  int rank_{};
  int size_{};

  size_t local_rows_{};
  size_t start_row_{};

  Image local_image_;
  std::vector<uint8_t> gathered_result_;
  Image result_image_;
};

}  // namespace morozov_n_sobels_filter
