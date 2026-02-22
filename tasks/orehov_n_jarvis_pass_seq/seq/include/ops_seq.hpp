#pragma once

#include "orehov_n_jarvis_pass_seq/common/include/common.hpp"
#include "task/include/task.hpp"

#include <vector>

namespace orehov_n_jarvis_pass_seq {

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
  
  [[nodiscard]] static double CheckLeft(Point a, Point b, Point c);
  [[nodiscard]] Point FindFirstElem() const;
  [[nodiscard]] static double Distance(Point a, Point b);
  [[nodiscard]] Point FindNext(Point current) const;

  std::vector<Point> res_;
  std::vector<Point> input_;
};

}  // namespace orehov_n_jarvis_pass_seq
