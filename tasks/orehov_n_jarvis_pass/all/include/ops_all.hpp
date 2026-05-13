#pragma once

#include <cstddef>
#include <vector>

#include "orehov_n_jarvis_pass/common/include/common.hpp"
#include "task/include/task.hpp"

namespace orehov_n_jarvis_pass {

class OrehovNJarvisPassALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit OrehovNJarvisPassALL(const InType &in);

  [[nodiscard]] static double CheckLeft(Point a, Point b, Point c);
  [[nodiscard]] static double Distance(Point a, Point b);

 private:
  struct BestState {
    Point point;
    bool valid = false;
  };

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] Point FindFirstElem() const;
  [[nodiscard]] Point FindNext(Point current) const;

  [[nodiscard]] static bool IsBetterPoint(const Point &current, const Point &candidate, const Point &best);
  [[nodiscard]] static BestState ReduceBestStates(const BestState &a, const BestState &b, const Point &current);

  [[nodiscard]] BestState LocalFindBest(const Point &current, size_t start, size_t end) const;
  [[nodiscard]] static BestState GlobalReduce(const std::vector<double> &all_data, int size, const Point &current);
  [[nodiscard]] static BestState FinalizeBestPoint(const double *global_data);

  std::vector<Point> res_;
  std::vector<Point> input_;
};

}  // namespace orehov_n_jarvis_pass
