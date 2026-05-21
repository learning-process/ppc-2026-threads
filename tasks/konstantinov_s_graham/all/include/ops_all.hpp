#pragma once
#include <cstddef>
#include <utility>
#include <vector>

#include "konstantinov_s_graham/common/include/common.hpp"
#include "task/include/task.hpp"

namespace konstantinov_s_graham {

class KonstantinovAGrahamALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit KonstantinovAGrahamALL(const InType &in);
  static constexpr double kKEps = 1e-10;

 private:
  int proc_rank_{0};
  int proc_num_{1};
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  static void RemoveDuplicates(std::vector<double> &xs, std::vector<double> &ys);

  [[nodiscard]] static bool IsLowerAnchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t lhs,
                                          size_t rhs);
  [[nodiscard]] static size_t FindAnchorIndex(const std::vector<double> &xs, const std::vector<double> &ys);

  [[nodiscard]] static double Dist2(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j);
  [[nodiscard]] static double CrossVal(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j,
                                       size_t k);

  static void FillIndexRange(std::vector<size_t> &idxs, size_t begin, size_t end, size_t anchor_idx);
  static void FillIndicesParallel(std::vector<size_t> &idxs, size_t point_count, size_t anchor_idx);

  [[nodiscard]] static std::vector<size_t> CollectAndSortIndices(const std::vector<double> &xs,
                                                                 const std::vector<double> &ys, size_t anchor_idx);

  [[nodiscard]] static bool CheckCollinearRange(const std::vector<double> &xs, const std::vector<double> &ys,
                                                size_t anchor_idx, const std::vector<size_t> &sorted_idxs, size_t begin,
                                                size_t end);

  [[nodiscard]] static bool AllCollinearWithAnchor(const std::vector<double> &xs, const std::vector<double> &ys,
                                                   size_t anchor_idx, const std::vector<size_t> &sorted_idxs);

  [[nodiscard]] static std::vector<std::pair<double, double>> BuildHullFromSorted(
      const std::vector<double> &xs, const std::vector<double> &ys, size_t anchor_idx,
      const std::vector<size_t> &sorted_idxs);

  [[nodiscard]] static std::vector<std::pair<double, double>> BuildHullFromCoords(const std::vector<double> &xs,
                                                                                  const std::vector<double> &ys);

  void ScatterInput(size_t total_size, std::vector<double> &local_xs, std::vector<double> &local_ys);
  void GatherLocalHull(const std::vector<std::pair<double, double>> &local_hull, std::vector<double> &gathered_xs,
                       std::vector<double> &gathered_ys) const;

  void BroadcastOutput();
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace konstantinov_s_graham
