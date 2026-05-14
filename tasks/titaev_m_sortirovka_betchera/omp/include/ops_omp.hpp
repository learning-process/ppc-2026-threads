#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevSortirovkaBetcheraOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit TitaevSortirovkaBetcheraOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static uint64_t PackDouble(double v) noexcept;
  static double UnpackDouble(uint64_t k) noexcept;

  static void LSDRadixSort(std::vector<double> &arr);
  static void BatcherOddEvenMerge(std::vector<double> &arr, size_t n);
  static inline void CompareSwap(std::vector<double> &arr, size_t i, size_t j) {
    if (arr[i] > arr[j]) {
      std::swap(arr[i], arr[j]);
    }
  }
};

}  // namespace titaev_m_sortirovka_betchera

#include "../src/ops_omp.cpp"
