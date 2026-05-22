#pragma once

#include <cstddef>

#include "lopatin_a_sobel_operator/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lopatin_a_sobel_operator {

class LopatinASobelOperatorSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit LopatinASobelOperatorSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::size_t h_ = 0;  // height
  std::size_t w_ = 0;  // width
};

}  // namespace lopatin_a_sobel_operator
