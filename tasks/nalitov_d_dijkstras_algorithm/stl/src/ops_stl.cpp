#include "nalitov_d_dijkstras_algorithm/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <thread>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "util/include/util.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

bool AddChecked(std::int64_t acc, Cost x, std::int64_t &out) {
  const auto xi = static_cast<std::int64_t>(x);
  if (xi > 0 && acc > std::numeric_limits<std::int64_t>::max() - xi) {
    return false;
  }
  if (xi < 0 && acc < std::numeric_limits<std::int64_t>::min() - xi) {
    return false;
  }
  out = acc + xi;
  return true;
}

bool FoldFinite(const std::vector<Cost> &row, OutType &sum) {
  std::int64_t acc = 0;
  for (Cost d : row) {
    if (d == kInf) {
      continue;
    }
    if (!AddChecked(acc, d, acc)) {
      return false;
    }
  }
  if (acc < 0) {
    return false;
  }
  sum = static_cast<OutType>(acc);
  return true;
}

bool ReduceShardResults(const std::vector<ShardResult> &shards, int worker_count, Cost *pick_cost,
                        NodeId *pick_vertex) {
  Cost pick_c = kInf;
  NodeId pick_v = -1;
  for (int shard_idx = 0; shard_idx < worker_count; ++shard_idx) {
    const ShardResult &sr = shards[static_cast<std::size_t>(shard_idx)];
    if (sr.id == -1) {
      continue;
    }
    const bool better_cost = sr.cost < pick_c;
    const bool tie_smaller_id = sr.cost == pick_c && (pick_v == -1 || sr.id < pick_v);
    if (better_cost || tie_smaller_id) {
      pick_c = sr.cost;
      pick_v = sr.id;
    }
  }
  *pick_cost = pick_c;
  *pick_vertex = pick_v;
  return pick_v != -1 && pick_c != kInf;
}

}  // namespace

NalitovDDijkstrasAlgorithmSTL::NalitovDDijkstrasAlgorithmSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

NalitovDDijkstrasAlgorithmSTL::~NalitovDDijkstrasAlgorithmSTL() {
  StopWorkers();
}

void NalitovDDijkstrasAlgorithmSTL::StopWorkers() {
  if (!pool_active_) {
    return;
  }
  {
    std::scoped_lock<std::mutex> hold(sync_);
    mode_ = WorkMode::kShutdown;
    ++epoch_;
    go_.notify_all();
  }
  for (std::thread &t : pool_) {
    if (t.joinable()) {
      t.join();
    }
  }
  pool_.clear();
  pool_active_ = false;
}

void NalitovDDijkstrasAlgorithmSTL::PartitionScanBest(int slot) {
  const int n = GetInput().n;
  Cost best_c = kInf;
  NodeId best_i = -1;
  for (int vid = slot; vid < n; vid += worker_slots_) {
    if (visited_[static_cast<std::size_t>(vid)] != 0) {
      continue;
    }
    const Cost d = dist_[static_cast<std::size_t>(vid)];
    if (d < best_c || (d == best_c && (best_i == -1 || vid < best_i))) {
      best_c = d;
      best_i = static_cast<NodeId>(vid);
    }
  }
  shard_results_[static_cast<std::size_t>(slot)] = ShardResult{.cost = best_c, .id = best_i};
}

void NalitovDDijkstrasAlgorithmSTL::PartitionPushDist(int slot) {
  const auto hub = static_cast<std::size_t>(pivot_);
  const auto &bundle = graph_[hub];
  const std::size_t m = bundle.size();
  if (m == 0) {
    return;
  }
  const std::size_t lo = (static_cast<std::size_t>(slot) * m) / static_cast<std::size_t>(worker_slots_);
  const std::size_t hi = ((static_cast<std::size_t>(slot) + 1U) * m) / static_cast<std::size_t>(worker_slots_);

  const Cost anchor = dist_[hub];
  if (anchor == kInf) {
    return;
  }

  for (std::size_t j = lo; j < hi; ++j) {
    const NodeId tgt = bundle[j].first;
    const Cost w = bundle[j].second;
    if (visited_[static_cast<std::size_t>(tgt)] != 0) {
      continue;
    }
    if (anchor > kInf - w) {
      continue;
    }
    const Cost cand = anchor + w;
    const std::size_t stripe = static_cast<std::size_t>(tgt) % kDistLockStripes;
    std::scoped_lock<std::mutex> hold(dist_stripes_.at(stripe));
    dist_[static_cast<std::size_t>(tgt)] = std::min(dist_[static_cast<std::size_t>(tgt)], cand);
  }
}

void NalitovDDijkstrasAlgorithmSTL::WorkerBody(int slot) {
  int seen = 0;
  for (;;) {
    WorkMode job = WorkMode::kParked;
    {
      std::unique_lock<std::mutex> hold(sync_);
      go_.wait(hold, [&] { return epoch_ > seen; });
      seen = epoch_;
      job = mode_;
    }
    if (job == WorkMode::kShutdown) {
      return;
    }
    if (job == WorkMode::kScanBest) {
      PartitionScanBest(slot);
    } else if (job == WorkMode::kPushDist) {
      PartitionPushDist(slot);
    }
    {
      std::scoped_lock<std::mutex> hold(sync_);
      if (--unfinished_ == 0) {
        settled_.notify_one();
      }
    }
  }
}

void NalitovDDijkstrasAlgorithmSTL::LaunchBarrier(WorkMode mode) {
  {
    std::scoped_lock<std::mutex> hold(sync_);
    mode_ = mode;
    unfinished_ = worker_slots_;
    ++epoch_;
    go_.notify_all();
  }
  std::unique_lock<std::mutex> hold(sync_);
  settled_.wait(hold, [&] { return unfinished_ == 0; });
}

bool NalitovDDijkstrasAlgorithmSTL::ValidationImpl() {
  if (GetOutput() != 0) {
    return false;
  }
  const InType &in = GetInput();
  constexpr int kCap = 10000;
  if (in.n <= 0 || in.n > kCap) {
    return false;
  }
  if (in.source < 0 || in.source >= in.n) {
    return false;
  }
  const auto ok = [&in](const Arc &a) {
    return a.from >= 0 && a.to >= 0 && a.from < in.n && a.to < in.n && a.weight >= 0;
  };
  return std::ranges::all_of(in.arcs, ok);
}

bool NalitovDDijkstrasAlgorithmSTL::PreProcessingImpl() {
  const InType &in = GetInput();
  graph_.assign(static_cast<std::size_t>(in.n), {});
  for (const Arc &arc : in.arcs) {
    graph_[static_cast<std::size_t>(arc.from)].emplace_back(arc.to, arc.weight);
  }

  dist_.assign(static_cast<std::size_t>(in.n), kInf);
  visited_.assign(static_cast<std::size_t>(in.n), 0);
  dist_[static_cast<std::size_t>(in.source)] = 0;

  worker_slots_ = std::max(1, ppc::util::GetNumThreads());
  shard_results_.assign(static_cast<std::size_t>(worker_slots_), ShardResult{});

  epoch_ = 0;
  unfinished_ = 0;
  mode_ = WorkMode::kParked;

  pool_.resize(static_cast<std::size_t>(worker_slots_));
  for (int thread_idx = 0; thread_idx < worker_slots_; ++thread_idx) {
    pool_[static_cast<std::size_t>(thread_idx)] =
        std::thread(&NalitovDDijkstrasAlgorithmSTL::WorkerBody, this, thread_idx);
  }
  pool_active_ = true;
  GetOutput() = 0;
  return true;
}

bool NalitovDDijkstrasAlgorithmSTL::RunImpl() {
  const InType &in = GetInput();
  if (graph_.size() != static_cast<std::size_t>(in.n)) {
    return false;
  }
  if (!pool_active_ || worker_slots_ <= 0) {
    return false;
  }

  for (int step = 0; step < in.n; ++step) {
    LaunchBarrier(WorkMode::kScanBest);
    Cost pick_c = kInf;
    NodeId pick_v = -1;
    if (!ReduceShardResults(shard_results_, worker_slots_, &pick_c, &pick_v)) {
      break;
    }
    visited_[static_cast<std::size_t>(pick_v)] = 1;
    pivot_ = pick_v;

    LaunchBarrier(WorkMode::kPushDist);
  }

  OutType total = 0;
  if (!FoldFinite(dist_, total)) {
    return false;
  }
  GetOutput() = total;
  return true;
}

bool NalitovDDijkstrasAlgorithmSTL::PostProcessingImpl() {
  StopWorkers();
  return GetOutput() >= 0;
}

}  // namespace nalitov_d_dijkstras_algorithm
