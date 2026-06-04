#pragma once

#include <cstdint>
#include <vector>

#include "egashin_k_radix_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace egashin_k_radix_simple_merge {

class EgashinKRadixSimpleMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit EgashinKRadixSimpleMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void CountingPass(const std::vector<uint64_t> &source, std::vector<uint64_t> &destination, int byte_index);
  static void RadixSort(std::vector<double> &data);

  OutType result_;
};

}  // namespace egashin_k_radix_simple_merge
