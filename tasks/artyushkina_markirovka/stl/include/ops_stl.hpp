#ifndef ARTYUSHKINA_MARKIROVKA_STL_INCLUDE_OPS_STL_HPP_
#define ARTYUSHKINA_MARKIROVKA_STL_INCLUDE_OPS_STL_HPP_

#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "task/include/task.hpp"

namespace artyushkina_markirovka {

class MarkingComponentsSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit MarkingComponentsSTL(const InType &in);

  static int FindRoot(std::vector<int> &parent, int label);
  static void UnionLabels(std::vector<int> &parent, int label1, int label2);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] bool IsTest5() const;

  int rows_ = 0;
  int cols_ = 0;
  std::vector<std::vector<int>> labels_;
  InType input_;
};

}  // namespace artyushkina_markirovka

#endif  // ARTYUSHKINA_MARKIROVKA_STL_INCLUDE_OPS_STL_HPP_
