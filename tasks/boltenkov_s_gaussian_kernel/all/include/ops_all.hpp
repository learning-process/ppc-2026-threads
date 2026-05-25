#pragma once

#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace boltenkov_s_gaussian_kernel {

class BoltenkovSGaussianKernelALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit BoltenkovSGaussianKernelALL(const InType &in);

 private:
  std::vector<std::vector<int>> kernel_;
  int shift_ = 4;

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void BcastSizes(int &n, int &m, int rank);
  static void ScatterRows(std::vector<std::vector<int>> &global_data, std::vector<std::vector<int>> &local_halo,
                          int &local_start_row, int &local_rows, int m, int rank, int size);
  void GatherResults(std::vector<std::vector<int>> &local_res, int local_start_row, int local_rows, int m, int rank,
                     int size);
  void GatherResultsRoot(std::vector<std::vector<int>> &local_res, int local_start_row, int local_rows, int m,
                         int size);
  static void GatherResultsOthers(std::vector<std::vector<int>> &local_res, int local_start_row, int local_rows, int m);

  std::vector<std::vector<int>> ApplyGaussianFilter(const std::vector<std::vector<int>> &local_halo,
                                                    int local_start_row, int local_rows, int m);

  static std::vector<std::vector<int>> CreateValidMatrix(int rows, int cols);
  static void FillTmpFromHalo(std::vector<std::vector<int>> &tmp, const std::vector<std::vector<int>> &local_halo,
                              int local_start_row, int local_rows, int m);

  static void SendRowsToOneProcess(int proc, int rows_per_proc, int n, int m,
                                   const std::vector<std::vector<int>> &global_data);

  static void FillLocalHaloForRoot(int local_start_row, int local_end_row, int n, int m,
                                   const std::vector<std::vector<int>> &global_data,
                                   std::vector<std::vector<int>> &local_halo);

  static void ReceiveRowsOnWorker(int m, int &local_start_row, int &local_rows,
                                  std::vector<std::vector<int>> &local_halo);
};

}  // namespace boltenkov_s_gaussian_kernel
