#pragma once

#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "votincev_d_radixmerge_sort/common/include/common.hpp"

namespace votincev_d_radixmerge_sort {

class VotincevDRadixMergeSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit VotincevDRadixMergeSortALL(InType in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  // ==============================
  // мои дополнительные функции ===
  static void LocalRadixSort(uint32_t *begin, uint32_t *end);
  static void Merge(const uint32_t *src, uint32_t *dst, int32_t left, int32_t mid, int32_t right);
  static void OmpLocalSortAndMerge(std::vector<uint32_t> &local_data);

  int32_t ScatterData(int32_t rank, int32_t n, int32_t local_n, const std::vector<int32_t> &send_counts,
                      const std::vector<int32_t> &displacements, std::vector<uint32_t> &local_data);

  void FinalMergeAndFormat(int32_t rank, int32_t size, int32_t n, int32_t min_val, std::vector<uint32_t> &gathered_data,
                           const std::vector<int32_t> &displacements);

  std::vector<int32_t> input_;
  std::vector<int32_t> output_;
};
}  // namespace votincev_d_radixmerge_sort
