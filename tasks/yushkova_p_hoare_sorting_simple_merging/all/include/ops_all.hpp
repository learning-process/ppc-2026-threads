#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

class YushkovaPHoareSortingSimpleMergingALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit YushkovaPHoareSortingSimpleMergingALL(const InType &in);

 private:
  static int HoarePartition(std::vector<int> &values, int left, int right);
  static void HoareQuickSort(std::vector<int> &values, int left, int right);
  static void SimpleMerge(const std::vector<int> &source, std::vector<int> &destination, std::size_t left,
                          std::size_t middle, std::size_t right);
  static void SortLocalStlParallel(std::vector<int> &values);
  static void MergeGatheredChunks(std::vector<int> &values, const std::vector<std::size_t> &chunk_sizes,
                                  const std::vector<std::size_t> &offsets);
  static void BroadcastVector(std::vector<int> &values, int rank);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace yushkova_p_hoare_sorting_simple_merging
