#pragma once

#include <queue>
#include <utility>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "task/include/task.hpp"

namespace artyushkina_markirovka {

struct NeighborOffset;

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

  static int FindRoot(int label);
  void UnionLabels(int label1, int label2);
  void BFS(int start_i, int start_j, int label);

  [[nodiscard]] bool IsValidNeighbor(int i, int j, const NeighborOffset &offset) const;
  void ProcessNeighbor(int i, int j, const NeighborOffset &offset, int label, std::queue<std::pair<int, int>> &q);

  int rows_ = 0;
  int cols_ = 0;
  std::vector<std::vector<int>> labels_;
  std::vector<int> equivalent_labels_;
};

}  // namespace artyushkina_markirovka
