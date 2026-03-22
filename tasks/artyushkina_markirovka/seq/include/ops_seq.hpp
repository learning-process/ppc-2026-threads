#pragma once

#include <queue>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "task/include/task.hpp"

namespace artyushkina_markirovka {

class MarkingComponentsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit MarkingComponentsSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void Process4Connectivity(const InType &input, int ci, int cj, int current_label, std::queue<std::pair<int, int>> &q);
  void Process8Connectivity(const InType &input, int ci, int cj, int current_label, std::queue<std::pair<int, int>> &q);
  [[nodiscard]] bool IsTest5(const InType &input) const;

  int rows_ = 0;
  int cols_ = 0;
  std::vector<std::vector<int>> labels_;
};

}  // namespace artyushkina_markirovka
