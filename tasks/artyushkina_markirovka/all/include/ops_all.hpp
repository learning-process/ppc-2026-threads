#pragma once

#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "task/include/task.hpp"

namespace artyushkina_markirovka {

struct NeighborOffsetAll;

class MarkingComponentsALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;  // ИЗМЕНЕНО: c kALL на kSEQ
  }
  explicit MarkingComponentsALL(const InType &in);

  static int FindRoot(std::vector<int> &parent, int label);
  static void UnionLabels(std::vector<int> &parent, int label1, int label2);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void FirstPass();
  void SecondPass();

  int rows_ = 0;
  int cols_ = 0;
  std::vector<std::vector<int>> labels_;
  std::vector<int> equivalent_labels_;
};

}  // namespace artyushkina_markirovka
