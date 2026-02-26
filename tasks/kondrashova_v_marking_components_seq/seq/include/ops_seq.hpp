#pragma once

#include <cstdint>
#include <vector>

#include "kondrashova_v_marking_components_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kondrashova_v_marking_components_seq {

class KondrashovaVTaskSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit KondrashovaVTaskSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int width_{};
  int height_{};
  std::vector<uint8_t> image_;
  std::vector<int> labels_1d_;
};

}  // namespace kondrashova_v_marking_components_seq