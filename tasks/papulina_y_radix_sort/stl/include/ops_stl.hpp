#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace papulina_y_radix_sort {

class PapulinaYRadixSortSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit PapulinaYRadixSortSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static uint64_t InBytes(double d);
  static double FromBytes(uint64_t bits);
  static void RadixSortParallel(double *arr, size_t size);
  static void ExecuteRadixPass(uint64_t *src, uint64_t *dst, size_t size, int byte_idx, unsigned int num_threads);
  static void CountFrequencies(const uint64_t *src, size_t size, int byte_idx, unsigned int num_threads,
                               std::vector<std::vector<size_t>> &local_hists);

  static void ComputeOffsets(unsigned int num_threads, const std::vector<std::vector<size_t>> &local_hists,
                             std::vector<std::vector<size_t>> &thread_pos);

  static void ReorderElements(const uint64_t *src, uint64_t *dst, size_t size, int byte_idx, unsigned int num_threads,
                              std::vector<std::vector<size_t>> &thread_pos);

  static std::pair<size_t, size_t> GetThreadRange(size_t size, unsigned int num_threads, unsigned int thread_idx);
  static const uint64_t kMask = 0x8000000000000000ULL;
};

}  // namespace papulina_y_radix_sort
