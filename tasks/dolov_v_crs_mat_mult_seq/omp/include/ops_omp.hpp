#pragma once

#include <vector>

#include "dolov_v_crs_mat_mult_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dolov_v_crs_mat_mult_seq {

struct RowData {
  std::vector<int> cols;
  std::vector<double> vals;
};

class DolovVCrsMatMultOmp : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit DolovVCrsMatMultOmp(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace dolov_v_crs_mat_mult_seq
