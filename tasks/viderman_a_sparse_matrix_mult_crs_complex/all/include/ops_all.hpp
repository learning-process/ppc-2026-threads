#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

class VidermanASparseMatrixMultCRSComplexALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit VidermanASparseMatrixMultCRSComplexALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void MasterScatter();
  void WorkerScatter();
  void ScatterData();

  void MasterGather();
  void WorkerGather();
  void GatherData();

  void SendDataProcesses();
  void ReceiveDataProcesses(CRSMatrix &out) const;
  void BroadcastData();

  void ComputeMultiply(int i, std::vector<int> &cols, std::vector<Complex> &vals);
  void MultiplyLocal();

  CRSMatrix local_a_;
  CRSMatrix b_;
  CRSMatrix local_c_;

  int a_rows_{0};
  int a_cols_{0};

  int rank_{0};
  int world_size_{1};
  int initialized_pipiline_{0};
};

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
