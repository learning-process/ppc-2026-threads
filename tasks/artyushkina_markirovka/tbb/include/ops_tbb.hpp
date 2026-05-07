#pragma once

#include <atomic>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "task/include/task.hpp"

namespace artyushkina_markirovka {

class MarkingComponentsTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit MarkingComponentsTBB(const InType &in);

  static int FindRoot(std::vector<int> &parent, int label);
  static void UnionLabels(std::vector<int> &parent, int label1, int label2);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] bool IsTest5() const;

  void ProcessFirstPass();
  void ResolveEquivalences();
  void RemapLabels();

  int rows_ = 0;
  int cols_ = 0;
  std::vector<std::vector<int>> labels_;
  std::vector<std::vector<int>> temp_labels_;
  std::vector<int> parent_;
  InType input_;
  std::atomic<int> next_label_{1};
  bool is_test5_ = false;
};

}  // namespace artyushkina_markirovka
