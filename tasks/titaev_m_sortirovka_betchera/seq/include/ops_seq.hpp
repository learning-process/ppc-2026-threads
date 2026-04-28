#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevSortirovkaBetcheraSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TitaevSortirovkaBetcheraSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ConvertToKeys(const InType &input, std::vector<uint64_t> &keys);
  static void RadixSort(std::vector<uint64_t> &keys);
  static void ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output);
  void BatcherSort();

  static void BatcherStep(OutType &result, size_t n, size_t step, size_t stage);
};

}  // namespace titaev_m_sortirovka_betchera
