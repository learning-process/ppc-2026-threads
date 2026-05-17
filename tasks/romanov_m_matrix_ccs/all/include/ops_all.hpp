#pragma once

#include <cstddef>
#include <vector>

#include "romanov_m_matrix_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace romanov_m_matrix_ccs {

class RomanovMMatrixCCSALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit RomanovMMatrixCCSALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void MultiplyColumn(size_t col_index, const MatrixCCS &a, const MatrixCCS &b, std::vector<double> &temp_v,
                             std::vector<size_t> &temp_r);

  static void SyncMatrixData(int rank, MatrixCCS &a, MatrixCCS &b);

  static void MasterCollect(int size, int chunk, int remainder, std::vector<std::vector<double>> &all_v,
                            std::vector<std::vector<size_t>> &all_r);

  static void WorkerSend(int local_count, std::vector<std::vector<double>> &local_v,
                         std::vector<std::vector<size_t>> &local_r);
};

}  // namespace romanov_m_matrix_ccs
