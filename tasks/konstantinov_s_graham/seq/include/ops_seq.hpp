#pragma once

#include "konstantinov_s_graham/common/include/common.hpp"
#include "task/include/task.hpp"

namespace konstantinov_s_graham {

class KonstantinovAGrahamSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KonstantinovAGrahamSEQ(const InType &in);
  static constexpr double kKEps = 1e-10;

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  static void RemoveDuplicates(std::vector<double> &xs, std::vector<double> &ys);
  static size_t FindAnchorIndex(const std::vector<double> &xs, const std::vector<double> &ys);
  static double Dist2(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j);
  double CrossVal(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j, size_t k);
  std::vector<size_t> CollectAndSortIndices(const std::vector<double> &xs, const std::vector<double> &ys,
                                            size_t anchor_idx);
  bool AllCollinearWithAnchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t anchor_idx,
                              const std::vector<size_t> &sorted_idxs);
  std::vector<std::pair<double, double>> BuildHullFromSorted(const std::vector<double> &xs,
                                                             const std::vector<double> &ys, size_t anchor_idx,
                                                             const std::vector<size_t> &sorted_idxs);
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace konstantinov_s_graham
