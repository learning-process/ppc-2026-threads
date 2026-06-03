#pragma once

#include <vector>

#include "krykov_e_sobel_op/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krykov_e_sobel_op {

class KrykovESobelOpSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit KrykovESobelOpSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> grayscale_;
  int width_ = 0;
  int height_ = 0;
};

}  // namespace krykov_e_sobel_op
