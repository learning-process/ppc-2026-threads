#pragma once

#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "task/include/task.hpp"

namespace artyushkina_markirovka {

struct NeighborOffsetAll;

class MarkingComponentsALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit MarkingComponentsALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static int FindRoot(std::vector<int> &parent, int label);
  static void UnionLabels(std::vector<int> &parent, int label1, int label2);
  
  void FirstPass();
  void SecondPass();

  [[nodiscard]] bool IsValidNeighbor(int i, int j, const NeighborOffsetAll &offset) const;
  void ProcessNeighborFirstPass(int i, int j, const NeighborOffsetAll &offset, 
                                std::vector<int> &neighbor_labels, int &min_label);

  int rows_ = 0;
  int cols_ = 0;
  std::vector<std::vector<int>> labels_;
  std::vector<int> equivalent_labels_;
};

}  // namespace artyushkina_markirovka