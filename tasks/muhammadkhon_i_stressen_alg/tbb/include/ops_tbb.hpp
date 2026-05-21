#pragma once

#include <cstddef>
#include <vector>

#include "muhammadkhon_i_stressen_alg/common/include/common.hpp"
#include "task/include/task.hpp"

namespace muhammadkhon_i_stressen_alg {

class MuhammadkhonIStressenAlgTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit MuhammadkhonIStressenAlgTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  size_t a_rows_ = 0;
  size_t a_cols_b_rows_ = 0;
  size_t b_cols_ = 0;
  size_t padded_n_ = 0;

  std::vector<double> padded_a_;
  std::vector<double> padded_b_;
  std::vector<double> result_c_;
};

}  // namespace muhammadkhon_i_stressen_alg
