#pragma once

#include <vector>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lukin_i_ench_contr_lin_hist {

class LukinITestTaskSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit LukinITestTaskSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void GetLocMinMax(std::vector<unsigned char> &loc_mins, std::vector<unsigned char> &loc_maxs);
  void Process(int min, int max);

  int thread_count_ = 0;
  int chunk_size_ = 0;
  int size_ = 0;
};

}  // namespace lukin_i_ench_contr_lin_hist
