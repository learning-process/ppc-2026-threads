#pragma once
#include <cstddef>
#include <utility>
#include <vector>

#include "konstantinov_s_graham/common/include/common.hpp"
#include "task/include/task.hpp"

namespace konstantinov_s_graham {

class KonstantinovAGrahamSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit KonstantinovAGrahamSTL(const InType &in);
  static constexpr double kKEps = 1e-10;

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  static bool IsLowerAnchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t lhs, size_t rhs);
  static size_t GetThreadCount(size_t n);
  static size_t FindLocalAnchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t begin, size_t end);
  static void RemoveDuplicates(std::vector<double> &xs, std::vector<double> &ys);
  static size_t FindAnchorIndex(const std::vector<double> &xs, const std::vector<double> &ys);
  static double Dist2(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j);
  static double CrossVal(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j, size_t k);
  static std::vector<size_t> CollectAndSortIndices(const std::vector<double> &xs, const std::vector<double> &ys,
                                                   size_t anchor_idx);
  static bool CheckCollinearRange(const std::vector<double> &xs, const std::vector<double> &ys, size_t anchor_idx,
                                  const std::vector<size_t> &sorted_idxs, size_t begin, size_t end);
  static void FillIndexRange(std::vector<size_t> &idxs, size_t begin, size_t end, size_t anchor_idx);
  static void FillIndicesParallel(std::vector<size_t> &idxs, size_t point_count, size_t anchor_idx,
                                  size_t thread_count);
  static bool AllCollinearWithAnchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t anchor_idx,
                                     const std::vector<size_t> &sorted_idxs);
  static std::vector<std::pair<double, double>> BuildHullFromSorted(const std::vector<double> &xs,
                                                                    const std::vector<double> &ys, size_t anchor_idx,
                                                                    const std::vector<size_t> &sorted_idxs);
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace konstantinov_s_graham
