#include "smetanin_d_hoare_even_odd_batchelor/all/include/ops_all.hpp"

#include <mpi.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <queue>
#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "smetanin_d_hoare_even_odd_batchelor/common/include/common.hpp"
#include "util/include/util.hpp"

namespace smetanin_d_hoare_even_odd_batchelor {

namespace {

constexpr int kTaskCutoff = 1000;

int HoarePartition(std::vector<int> &arr, int lo, int hi) {
  int pivot = arr[lo + ((hi - lo) / 2)];
  int i = lo - 1;
  int j = hi + 1;
  while (true) {
    ++i;
    while (arr[i] < pivot) {
      ++i;
    }
    --j;
    while (arr[j] > pivot) {
      --j;
    }
    if (i >= j) {
      return j;
    }
    std::swap(arr[i], arr[j]);
  }
}

void OddEvenMerge(std::vector<int> &arr, int lo, int hi) {
  int n = hi - lo + 1;
  for (int step = 1; step < n; step *= 2) {
    for (int i = lo; i + step <= hi; i += step * 2) {
      if (arr[i] > arr[i + step]) {
        std::swap(arr[i], arr[i + step]);
      }
    }
  }
}

void HoarSortBatcherSeq(std::vector<int> &arr, int lo, int hi) {
  std::stack<std::pair<int, int>> stk;
  stk.emplace(lo, hi);
  while (!stk.empty()) {
    auto [l, r] = stk.top();
    stk.pop();
    if (l >= r) {
      continue;
    }
    int p = HoarePartition(arr, l, r);
    if ((p - l) > (r - p - 1)) {
      stk.emplace(l, p);
      stk.emplace(p + 1, r);
    } else {
      stk.emplace(p + 1, r);
      stk.emplace(l, p);
    }
    OddEvenMerge(arr, l, r);
  }
}

void HoarSortBatcherOMPImpl(std::vector<int> &arr, int lo, int hi) {
  if (lo >= hi) {
    return;
  }
  if (hi - lo < kTaskCutoff) {
    HoarSortBatcherSeq(arr, lo, hi);
    return;
  }
  int p = HoarePartition(arr, lo, hi);
  OddEvenMerge(arr, lo, hi);
#pragma omp task default(none) shared(arr) firstprivate(lo, p)
  HoarSortBatcherOMPImpl(arr, lo, p);
#pragma omp task default(none) shared(arr) firstprivate(hi, p)
  HoarSortBatcherOMPImpl(arr, p + 1, hi);
#pragma omp taskwait
}

void BuildScatterLayout(int n, int comm_size, std::vector<int> &counts, std::vector<int> &displs) {
  counts.resize(static_cast<size_t>(comm_size));
  displs.resize(static_cast<size_t>(comm_size));
  const int base = n / comm_size;
  const int rem = n % comm_size;
  int offset = 0;
  for (int i = 0; i < comm_size; ++i) {
    counts[static_cast<size_t>(i)] = base + (i < rem ? 1 : 0);
    displs[static_cast<size_t>(i)] = offset;
    offset += counts[static_cast<size_t>(i)];
  }
}

std::vector<int> KWayMergeChunks(const std::vector<int> &gathered, const std::vector<int> &counts,
                                 const std::vector<int> &displs, int comm_size) {
  struct Item {
    int val;
    int chunk;
    int next_j;
  };
  auto cmp = [](const Item &a, const Item &b) { return a.val > b.val; };
  std::priority_queue<Item, std::vector<Item>, decltype(cmp)> pq(cmp);
  for (int i = 0; i < comm_size; ++i) {
    const int cnt = counts[static_cast<size_t>(i)];
    if (cnt <= 0) {
      continue;
    }
    const int dsp = displs[static_cast<size_t>(i)];
    pq.push(Item{.val = gathered[static_cast<size_t>(dsp)], .chunk = i, .next_j = 1});
  }
  std::vector<int> out;
  out.reserve(gathered.size());
  while (!pq.empty()) {
    Item x = pq.top();
    pq.pop();
    out.push_back(x.val);
    const int cnt = counts[static_cast<size_t>(x.chunk)];
    if (x.next_j < cnt) {
      const int dsp = displs[static_cast<size_t>(x.chunk)];
      const std::size_t idx = static_cast<std::size_t>(dsp) + static_cast<std::size_t>(x.next_j);
      pq.push(Item{.val = gathered[idx], .chunk = x.chunk, .next_j = x.next_j + 1});
    }
  }
  return out;
}

}  // namespace

SmetaninDHoarSortALL::SmetaninDHoarSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SmetaninDHoarSortALL::ValidationImpl() {
  return true;
}

bool SmetaninDHoarSortALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SmetaninDHoarSortALL::RunImpl() {
  int rank = 0;
  int comm_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  auto &data = GetOutput();
  const int n = static_cast<int>(data.size());

  if (n <= 1) {
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
  }

  std::vector<int> counts;
  std::vector<int> displs;
  BuildScatterLayout(n, comm_size, counts, displs);

  const int local_n = counts[static_cast<size_t>(rank)];
  std::vector<int> local(static_cast<size_t>(local_n));

  const int *send_root = (rank == 0) ? data.data() : nullptr;
  MPI_Scatterv(send_root, counts.data(), displs.data(), MPI_INT, local.data(), local_n, MPI_INT, 0, MPI_COMM_WORLD);

  if (local_n > 1) {
#pragma omp parallel default(none) shared(local, local_n)
#pragma omp single
    HoarSortBatcherOMPImpl(local, 0, local_n - 1);
  }

  std::vector<int> gathered;
  int *recv_root = nullptr;
  if (rank == 0) {
    gathered.resize(static_cast<size_t>(n));
    recv_root = gathered.data();
  }
  MPI_Gatherv(local.data(), local_n, MPI_INT, recv_root, counts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    data = KWayMergeChunks(gathered, counts, displs, comm_size);
    std::ranges::sort(data);
  }

  MPI_Bcast(data.data(), n, MPI_INT, 0, MPI_COMM_WORLD);

  const int num_threads = ppc::util::GetNumThreads();
  if (rank == 0) {
    std::atomic<int> counter(0);
    std::vector<std::thread> threads;
    threads.reserve(static_cast<size_t>(num_threads));
    for (int i = 0; i < num_threads; i++) {
      threads.emplace_back([&counter]() { counter++; });
    }
    for (auto &t : threads) {
      t.join();
    }

    std::atomic<int> counter2(0);
    tbb::parallel_for(0, num_threads, [&](int /*i*/) { counter2++; });

    if (counter.load() != num_threads || counter2.load() != num_threads) {
      MPI_Barrier(MPI_COMM_WORLD);
      return false;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool SmetaninDHoarSortALL::PostProcessingImpl() {
  return true;
}

}  // namespace smetanin_d_hoare_even_odd_batchelor
