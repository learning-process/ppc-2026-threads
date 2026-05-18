#pragma once

#include "sizov_d_sparse_crs_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sizov_d_sparse_crs_mult {

class SizovDSparseCRSMultOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit SizovDSparseCRSMultOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sizov_d_sparse_crs_mult
