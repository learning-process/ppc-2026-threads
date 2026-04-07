#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "Nazarova_K_rad_sort_batcher_metod/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nazarova_k_rad_sort_batcher_metod_processes {

class NazarovaKRadSortBatcherMetodSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit NazarovaKRadSortBatcherMetodSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void LSDRadixSort(std::vector<double> &array);
  static std::uint64_t PackDouble(double v) noexcept;
  static double UnpackDouble(std::uint64_t k) noexcept;
  static void BatcherOddEvenMerge(std::vector<double> &array, std::size_t n);
  static void BatcherMergeSort(std::vector<double> &array);
  static void BlocksComparing(std::vector<double> &arr, std::size_t i, std::size_t step);
};

}  // namespace nazarova_k_rad_sort_batcher_metod_processes
