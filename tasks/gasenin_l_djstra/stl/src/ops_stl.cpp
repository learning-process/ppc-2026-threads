#include "gasenin_l_djstra/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <thread>
#include <vector>

#include "gasenin_l_djstra/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gasenin_l_djstra {

GaseninLDjstraSTL::GaseninLDjstraSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool GaseninLDjstraSTL::ValidationImpl() {
  return GetInput() > 0;
}

bool GaseninLDjstraSTL::PreProcessingImpl() {
  const InType n = GetInput();
  const InType inf = std::numeric_limits<InType>::max();

  dist_.assign(n, inf);
  visited_.assign(n, 0);
  dist_[0] = 0;

  num_threads_ = ppc::util::GetNumThreads();
  local_min_.assign(num_threads_, inf);
  local_vert_.assign(num_threads_, -1);

  generation_ = 0;
  pending_ = 0;
  phase_ = Phase::kIdle;

  workers_.resize(num_threads_);
  for (int thread_id = 0; thread_id < num_threads_; ++thread_id) {
    workers_[thread_id] = std::thread(&GaseninLDjstraSTL::WorkerLoop, this, thread_id);
  }

  return true;
}

void GaseninLDjstraSTL::DoFindMin(int thread_id) {
  const InType n = GetInput();
  const InType inf = std::numeric_limits<InType>::max();

  InType thread_min = inf;
  InType thread_vert = -1;
  for (int idx = thread_id; idx < n; idx += num_threads_) {
    if (visited_[idx] == 0 && dist_[idx] < thread_min) {
      thread_min = dist_[idx];
      thread_vert = idx;
    }
  }
  local_min_[thread_id] = thread_min;
  local_vert_[thread_id] = thread_vert;
}

void GaseninLDjstraSTL::DoRelax(int thread_id) {
  const InType n = GetInput();
  const InType inf = std::numeric_limits<InType>::max();
  const InType src = global_vertex_;

  for (int vertex = thread_id; vertex < n; vertex += num_threads_) {
    if (visited_[vertex] == 0 && vertex != src && dist_[src] != inf) {
      dist_[vertex] = std::min(dist_[vertex], dist_[src] + std::abs(src - vertex));
    }
  }
}

void GaseninLDjstraSTL::WorkerLoop(int thread_id) {
  int my_gen = 0;

  while (true) {
    Phase current_phase{};
    {
      std::unique_lock<std::mutex> lock(mtx_);
      cv_start_.wait(lock, [&] { return generation_ > my_gen; });
      my_gen = generation_;
      current_phase = phase_;
    }

    if (current_phase == Phase::kStop) {
      return;
    }

    if (current_phase == Phase::kFindMin) {
      DoFindMin(thread_id);
    } else {
      DoRelax(thread_id);
    }

    {
      std::scoped_lock<std::mutex> lock(mtx_);
      if (--pending_ == 0) {
        cv_done_.notify_one();
      }
    }
  }
}

void GaseninLDjstraSTL::Dispatch(Phase phase) {
  {
    std::scoped_lock<std::mutex> lock(mtx_);
    phase_ = phase;
    pending_ = num_threads_;
    ++generation_;
    cv_start_.notify_all();
  }
  std::unique_lock<std::mutex> lock(mtx_);
  cv_done_.wait(lock, [&] { return pending_ == 0; });
}

bool GaseninLDjstraSTL::RunImpl() {
  const InType n = GetInput();
  const InType inf = std::numeric_limits<InType>::max();

  for (int iter = 0; iter < n; ++iter) {
    Dispatch(Phase::kFindMin);

    InType global_min = inf;
    global_vertex_ = -1;
    for (int thread_idx = 0; thread_idx < num_threads_; ++thread_idx) {
      if (local_min_[thread_idx] < global_min) {
        global_min = local_min_[thread_idx];
        global_vertex_ = local_vert_[thread_idx];
      }
    }

    if (global_vertex_ == -1 || global_min == inf) {
      break;
    }

    visited_[global_vertex_] = 1;
    std::ranges::fill(local_min_, inf);
    std::ranges::fill(local_vert_, -1);

    Dispatch(Phase::kRelax);
  }

  int64_t total_sum = 0;
  for (int idx = 0; idx < n; ++idx) {
    if (dist_[idx] != inf) {
      total_sum += dist_[idx];
    }
  }
  GetOutput() = static_cast<OutType>(total_sum);
  return true;
}

bool GaseninLDjstraSTL::PostProcessingImpl() {
  {
    std::scoped_lock<std::mutex> lock(mtx_);
    phase_ = Phase::kStop;
    ++generation_;
    cv_start_.notify_all();
  }
  for (auto &worker : workers_) {
    worker.join();
  }
  workers_.clear();
  return true;
}

}  // namespace gasenin_l_djstra
