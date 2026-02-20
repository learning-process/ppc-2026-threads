#pragma once

#include <vector>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

class DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  // Локальный вектор для сортировки
  std::vector<double> local_data_;
};

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge