#pragma once

#include <mpi.h>

#include <cstddef>
#include <utility>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "task/include/task.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

class PerepelkinIConvexHullGrahamScanALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit PerepelkinIConvexHullGrahamScanALL(const InType &in);

 private:
  int proc_rank_{};
  int proc_num_{};
  MPI_Datatype MPI_POINT_{};

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  // Processes parallelization
  void BcastSizes(size_t &original_size, size_t &padded_size);
  void DistributeData(const size_t &original_size, const size_t &padded_size,
                      std::vector<std::pair<double, double>> &local_data);
  void BcastOutput();
  void FindGlobalPivot(std::vector<std::pair<double, double>> &local_data, const std::pair<double, double> &local_pivot,
                       std::pair<double, double> &global_pivot);
  void GatherSortedBlocks(std::vector<std::pair<double, double>> &local_data, std::vector<int> &real_counts,
                          std::vector<int> &real_displs, std::vector<std::pair<double, double>> &gathered);

  // Merge blocks
  std::vector<std::pair<double, double>> MergeSortedBlocks(const std::vector<std::pair<double, double>> &gathered,
                                                           const std::vector<int> &counts,
                                                           const std::vector<int> &displs,
                                                           const std::pair<double, double> &pivot);
  std::vector<std::pair<double, double>> MergeSortedBlocksParallel(
      const std::vector<std::pair<double, double>> &gathered, const std::vector<int> &counts,
      const std::vector<int> &displs, const std::pair<double, double> &pivot);
  std::vector<std::pair<double, double>> MergeBlocksRange(
      const std::vector<std::vector<std::pair<double, double>>> &blocks, int left, int right,
      const std::pair<double, double> &pivot);
  static std::vector<std::pair<double, double>> MergeTwoBlocks(const std::vector<std::pair<double, double>> &left,
                                                               const std::vector<std::pair<double, double>> &right,
                                                               const std::pair<double, double> &pivot);

  // Threads parallelization
  static void FindPivotParallel(const std::vector<std::pair<double, double>> &pts,
                                std::pair<double, double> &local_pivot);
  static void ParallelSort(std::vector<std::pair<double, double>> &data, const std::pair<double, double> &pivot);

  // Sequential
  static void HullConstruction(std::vector<std::pair<double, double>> &hull,
                               const std::vector<std::pair<double, double>> &pts,
                               const std::pair<double, double> &pivot);

  // Helpers
  [[nodiscard]] static bool IsBetterPivot(const std::pair<double, double> &a, const std::pair<double, double> &b);
  static bool AngleCmp(const std::pair<double, double> &a, const std::pair<double, double> &b,
                       const std::pair<double, double> &pivot);
  static double Orientation(const std::pair<double, double> &p, const std::pair<double, double> &q,
                            const std::pair<double, double> &r);
};

}  // namespace perepelkin_i_convex_hull_graham_scan
