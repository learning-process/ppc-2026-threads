#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

class YushkovaPHoareSortingSimpleMergingTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit YushkovaPHoareSortingSimpleMergingTBB(const InType &in);

 private:
  static int HoarePartition(std::vector<int> &values, int left, int right);
  static void HoareQuickSort(std::vector<int> &values, int left, int right);
  static void SimpleMerge(const std::vector<int> &source, std::vector<int> &destination, size_t left, size_t middle,
                          size_t right);
  static void SortBlockIfNeeded(std::vector<int> &data, size_t size, size_t block_index);
  static void SortBlocks(std::vector<int> &data, size_t size);
  static void MergeChunk(const std::vector<int> &source, std::vector<int> &destination, size_t size, size_t merge_width,
                         size_t merge_index);
  static void MergePass(std::vector<int> &data, size_t size, size_t merge_width);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> data_;
};

}  // namespace yushkova_p_hoare_sorting_simple_merging
