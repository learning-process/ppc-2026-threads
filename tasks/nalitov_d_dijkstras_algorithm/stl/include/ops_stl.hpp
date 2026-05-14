#pragma once

#include <array>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nalitov_d_dijkstras_algorithm {

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4324)
#endif
struct alignas(64) ShardResult {
  Cost cost = kInf;
  NodeId id = -1;
};
#ifdef _MSC_VER
#  pragma warning(pop)
#endif

class NalitovDDijkstrasAlgorithmSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit NalitovDDijkstrasAlgorithmSTL(const InType &in);
  ~NalitovDDijkstrasAlgorithmSTL() override;

 private:
  static constexpr std::size_t kDistLockStripes = 256;

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void StopWorkers();

  enum class WorkMode : std::uint8_t {
    kParked,
    kScanBest,
    kPushDist,
    kShutdown,
  };

  void WorkerBody(int slot);
  void PartitionScanBest(int slot);
  void PartitionPushDist(int slot);
  void LaunchBarrier(WorkMode mode);

  using OutgoingTable = std::vector<std::vector<std::pair<NodeId, Cost>>>;
  OutgoingTable graph_;

  std::vector<Cost> dist_;
  std::vector<char> visited_;

  int worker_slots_{};
  std::vector<std::thread> pool_;
  std::vector<ShardResult> shard_results_;

  NodeId pivot_{};

  std::mutex sync_;
  std::array<std::mutex, kDistLockStripes> dist_stripes_{};
  std::condition_variable go_;
  std::condition_variable settled_;
  WorkMode mode_{WorkMode::kParked};
  int epoch_{};
  int unfinished_{};

  bool pool_active_{false};
};

}  // namespace nalitov_d_dijkstras_algorithm
