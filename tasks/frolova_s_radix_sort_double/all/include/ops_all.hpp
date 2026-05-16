#pragma once

#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace frolova_s_radix_sort_double {

class FrolovaSRadixSortDoubleALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit FrolovaSRadixSortDoubleALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixSortChunk(std::vector<double> &chunk);
  static void ProcessChunk(std::vector<double> &chunk);
  static std::vector<double> SimpleMerge(const std::vector<double> &a,
                                         const std::vector<double> &b);
  static std::vector<double> ParallelMerge(std::vector<std::vector<double>> &chunks);
};

}  // namespace frolova_s_radix_sort_double