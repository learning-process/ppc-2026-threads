#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevSortirovkaBetcheraSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit TitaevSortirovkaBetcheraSTL(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  static void ConvertToKeys(const InType &input, std::vector<uint64_t> &keys);

  static void RadixSortParallel(std::vector<uint64_t> &keys);

  static void ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output);

  void BatcherSortParallel();

  static void BatcherStepParallel(OutType &res, size_t n, size_t step, size_t stage);
};

}  // namespace titaev_m_sortirovka_betchera
