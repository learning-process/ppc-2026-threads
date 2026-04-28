#pragma once

#include <cstdint>
#include <vector>

#include "pikhotskiy_r_vertical_gauss_filter/common/include/common.hpp"
#include "task/include/task.hpp"

namespace pikhotskiy_r_vertical_gauss_filter {

class PikhotskiyRVerticalGaussFilterSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit PikhotskiyRVerticalGaussFilterSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  void RunVerticalPassForStripe(int x_begin, int x_end);
  void RunHorizontalPassForStripe(int x_begin, int x_end);

  int width_ = 0;
  int height_ = 0;
  int stripe_width_ = 1;
  std::vector<std::uint8_t> source_;
  std::vector<int> vertical_buffer_;
  std::vector<std::uint8_t> result_buffer_;
};

}  // namespace pikhotskiy_r_vertical_gauss_filter
