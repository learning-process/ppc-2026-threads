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

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ConvertToKeys(const InType &input, std::vector<uint64_t> &keys);
  static void RadixSort(std::vector<uint64_t> &keys);
  static void CountSequential(const std::vector<uint64_t> &keys, std::vector<std::size_t> &count, int pass);
  static void CountParallel(const std::vector<uint64_t> &keys, std::vector<std::size_t> &count, int pass,
                            unsigned int num_threads);
  static void RadixCountPass(std::vector<uint64_t> &keys, std::vector<uint64_t> &tmp, int pass);
  static void ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output);
  void BatcherSort();
  static void BatcherStage(OutType &result, std::size_t array_size, std::size_t block, std::size_t step);
};

}  // namespace titaev_m_sortirovka_betchera
