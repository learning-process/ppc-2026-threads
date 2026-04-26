#pragma once

#include <vector>

#include "orehov_n_jarvis_pass/common/include/common.hpp"
#include "task/include/task.hpp"

namespace orehov_n_jarvis_pass {

class OrehovNJarvisPassTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit OrehovNJarvisPassTBB(const InType &in);

  [[nodiscard]] static double CheckLeft(Point a, Point b, Point c);
  [[nodiscard]] static double Distance(Point a, Point b);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] Point FindFirstElem() const;
  [[nodiscard]] Point FindNext(Point current) const;

  std::vector<Point> res_;
  std::vector<Point> input_;
};

}  // namespace orehov_n_jarvis_pass
