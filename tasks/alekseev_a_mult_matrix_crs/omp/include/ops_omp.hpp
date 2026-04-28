#pragma once

#include <cstddef>
#include <vector>

#include "alekseev_a_mult_matrix_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace alekseev_a_mult_matrix_crs {

class AlekseevAMultMatrixCRSOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit AlekseevAMultMatrixCRSOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ProcessRow(std::size_t i, const CRSMatrix &a, const CRSMatrix &b, std::vector<double> &temp_v,
                         std::vector<std::size_t> &temp_c, std::vector<double> &accum, std::vector<int> &touched_flag,
                         std::vector<std::size_t> &touched_cols);
};

}  // namespace alekseev_a_mult_matrix_crs
