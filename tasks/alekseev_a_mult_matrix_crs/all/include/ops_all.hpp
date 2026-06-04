#pragma once

#include <cstddef>
#include <vector>

#include "alekseev_a_mult_matrix_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace alekseev_a_mult_matrix_crs {

class AlekseevAMultMatrixCRSALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit AlekseevAMultMatrixCRSALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void BroadcastData(CRSMatrix &a, CRSMatrix &b, int rank);

  static void GatherData(CRSMatrix &out, const InType &input, const std::vector<std::vector<double>> &local_v,
                         const std::vector<std::vector<std::size_t>> &local_c, const std::vector<int> &send_counts,
                         const std::vector<int> &displs, int rank, int size);
  // ...
  static void ProcessRow(std::size_t i, const CRSMatrix &a, const CRSMatrix &b, std::vector<double> &temp_v,
                         std::vector<std::size_t> &temp_c, std::vector<double> &accum, std::vector<int> &touched_flag,
                         std::vector<std::size_t> &touched_cols);
};

}  // namespace alekseev_a_mult_matrix_crs
