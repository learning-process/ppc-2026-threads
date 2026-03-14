#pragma once

#include <vector>

#include "orehov_n_jarvis_pass/common/include/common.hpp"
#include "task/include/task.hpp"

namespace orehov_n_jarvis_pass {

class OrehovNJarvisPassOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit OrehovNJarvisPassOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] static double CheckLeft(Point a, Point b, Point c);
  [[nodiscard]] Point FindFirstElem() const;
  [[nodiscard]] static double Distance(Point a, Point b);
  [[nodiscard]] Point FindNext(Point current) const;
  [[nodiscard]] Point FindLocalBest(Point current, Point initial_next) const;
  [[nodiscard]] void UpdateGlobalBest(Point current, Point local_next, Point &global_next) const;

  std::vector<Point> res_;
  std::vector<Point> input_;
};

}  // namespace orehov_n_jarvis_pass
