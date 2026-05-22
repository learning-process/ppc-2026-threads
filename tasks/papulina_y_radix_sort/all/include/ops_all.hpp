#pragma once

#include <cstdint>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace papulina_y_radix_sort {

class PapulinaYRadixSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit PapulinaYRadixSortALL(const InType &in);

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
  static std::vector<double> SimpleMerge(const std::vector<double> &a, const std::vector<double> &b);
};

}  // namespace papulina_y_radix_sort
