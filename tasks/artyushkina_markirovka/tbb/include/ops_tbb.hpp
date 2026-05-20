#ifndef ARTYUSHKINA_MARKIROVKA_TBB_INCLUDE_OPS_TBB_HPP_
#define ARTYUSHKINA_MARKIROVKA_TBB_INCLUDE_OPS_TBB_HPP_

#include <tbb/spin_mutex.h>

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

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int FindRoot(int label);
  void UnionLabels(int label1, int label2);

  void InitLabelsTbb();
  void MergeHorizontalPairsTbb();
  void MergeVerticalPairsTbb();
  void FinalizeRootsTbb();
  void NormalizeLabelsTbb();

  int rows_ = 0;
  int cols_ = 0;
  std::vector<int> labels_;
  std::vector<int> parent_;
  InType input_;
  int current_label_ = 0;

  mutable tbb::spin_mutex dsu_mutex_;
};

}  // namespace artyushkina_markirovka

#endif  // ARTYUSHKINA_MARKIROVKA_TBB_INCLUDE_OPS_TBB_HPP_
