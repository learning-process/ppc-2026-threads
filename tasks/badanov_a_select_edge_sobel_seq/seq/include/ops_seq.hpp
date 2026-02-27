#pragma once

#include <vector>
#include <cstdint>

#include "example_threads/common/include/common.hpp"
#include "task/include/task.hpp"

namespace badanov_a_select_edge_sobel {

class BadanovASelectEdgeSobelSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BadanovASelectEdgeSobelSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static constexpr std::array<std::array<int, 3>, 3> KERNEL_X = {{
    {{-1, 0, 1}},
    {{-2, 0, 2}},
    {{-1, 0, 1}}
  }};
  
  static constexpr std::array<std::array<int, 3>, 3> KERNEL_Y = {{
    {{-1, -2, -1}},
    {{0, 0, 0}},
    {{1, 2, 1}}
  }};

  int width_ = 0;
  int height_ = 0;
  int threshold_ = 50;
};

}