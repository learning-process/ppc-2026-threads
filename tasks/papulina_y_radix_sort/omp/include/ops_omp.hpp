#pragma once
#include <cstdint>
#include <span>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"
namespace papulina_y_radix_sort {

class PapulinaYRadixSortOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit PapulinaYRadixSortOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static const uint64_t kMask = 0x8000000000000000ULL;
  static void SortByByte(uint64_t *bytes, uint64_t *out, int byte, int size);
  static uint64_t InBytes(double d);
  static double FromBytes(uint64_t bits);
  static void RadixSort(double *arr, int size);
  void MergeChunks(std::vector<std::span<double>> &chunks, double *result);
  static std::vector<double> Merge(const std::vector<double> &a, const std::vector<double> &b);
};

}  // namespace papulina_y_radix_sort
