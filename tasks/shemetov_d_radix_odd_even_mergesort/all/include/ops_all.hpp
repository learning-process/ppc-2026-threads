#pragma once

#include <cstddef>
#include <vector>

#include "shemetov_d_radix_odd_even_mergesort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shemetov_d_radix_odd_even_mergesort {

class ShemetovDRadixOddEvenMergeSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ShemetovDRadixOddEvenMergeSortALL(const InType &in);

 private:
  static void ScatterData(std::vector<int> &global_array, std::vector<int> &local_array, size_t chunk, int rank,
                          int ranks, int mpi_init);
  static void MPISort(std::vector<int> &local_array, size_t chunk);
  void MPIMerge(size_t chunk);
  static void GatherData(std::vector<int> &global_array, std::vector<int> &local_array, size_t chunk, int rank,
                         int ranks, int mpi_init);

  static void RadixSort(std::vector<int> &array, size_t left, size_t right, std::vector<int> &buffer,
                        std::vector<int> &position);
  static void OddEvenMerge(std::vector<int> &array, size_t start, size_t segment);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> array_;

  int offset_{0};
  size_t size_{0};
  size_t power_{0};
};

}  // namespace shemetov_d_radix_odd_even_mergesort
