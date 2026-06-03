#pragma once

#include <vector>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace baldin_a_radix_sort {

class BaldinARadixSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit BaldinARadixSortALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void DistributeData(int rank, int size, int &n);
  void GatherData(int rank);

  void LocalSort(int num_threads);
  void GlobalMerge(int num_threads, int size, int n);

  std::vector<int> counts_;
  std::vector<int> displs_;
  std::vector<int> local_data_;
};

}  // namespace baldin_a_radix_sort
