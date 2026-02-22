#pragma once

#include "orehov_n_Jarvis_pass_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace orehov_n_Jarvis_pass_seq {

class OrehovNJarvisPassSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit OrehovNJarvisPassSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  
  double CheckLeft(Point A, Point B, Point C) const;
  Point FindFirstElem() const;
  double distance(Point A, Point B) const;

  std::vector<Point> res;
  std::vector<Point> input;
};

}  // namespace orehov_n_Jarvis_pass_seq
