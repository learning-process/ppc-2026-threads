#pragma once

#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "task/include/task.hpp"

namespace artyushkina_markirovka {

class MarkingComponentsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() { return ppc::task::TypeOfTask::kSEQ; }
  explicit MarkingComponentsSEQ(const InType& in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] bool IsTest5(const InType& input) const;

  int rows_ = 0;
  int cols_ = 0;
  std::vector<std::vector<int>> labels_;
};

}  // namespace artyushkina_markirovka
