#pragma once

#include <cstdint>
#include <vector>

#include "kondrashova_v_marking_components/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kondrashova_v_marking_components {

class KondrashovaVTaskALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KondrashovaVTaskALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int width_{};
  int height_{};
  int rank_{};
  int world_size_{};
  std::vector<uint8_t> image_;
  std::vector<int> labels_1d_;
};

}  // namespace kondrashova_v_marking_components
