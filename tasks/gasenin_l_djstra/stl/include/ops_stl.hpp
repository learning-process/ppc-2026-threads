#pragma once

#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>
#include <vector>

#include "gasenin_l_djstra/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gasenin_l_djstra {

class GaseninLDjstraSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit GaseninLDjstraSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  enum class Phase : std::uint8_t {
    kIdle,
    kFindMin,
    kRelax,
    kStop,
  };

  void WorkerLoop(int thread_id);
  void DoFindMin(int thread_id);
  void DoRelax(int thread_id);
  void Dispatch(Phase phase);

  std::vector<InType> dist_;
  std::vector<char> visited_;
  std::vector<InType> local_min_;
  std::vector<InType> local_vert_;
  InType global_vertex_ = -1;
  int num_threads_ = 0;

  std::mutex mtx_;
  std::condition_variable cv_start_;
  std::condition_variable cv_done_;
  Phase phase_ = Phase::kIdle;
  int generation_ = 0;
  int pending_ = 0;

  std::vector<std::thread> workers_;
};

}  // namespace gasenin_l_djstra
