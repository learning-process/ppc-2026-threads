#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevSortirovkaBetcheraTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit TitaevSortirovkaBetcheraTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ConvertToKeys(const InType &input, std::vector<uint64_t> &keys);
  static void RadixSort(std::vector<uint64_t> &keys);
  static void ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output);
  void BatcherSort();
  static void BatcherStage(OutType &result, std::size_t array_size, std::size_t block, std::size_t step);
};

}  // namespace titaev_m_sortirovka_betchera
