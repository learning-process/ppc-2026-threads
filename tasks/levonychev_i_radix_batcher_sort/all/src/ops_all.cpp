#include "levonychev_i_radix_batcher_sort/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <future>
#include <iterator>
#include <thread>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"

namespace levonychev_i_radix_batcher_sort {

LevonychevIRadixBatcherSortALL::LevonychevIRadixBatcherSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void LevonychevIRadixBatcherSortALL::LocalRadixSort(std::vector<int> &arr) {
  if (arr.empty()) {
    return;
  }
  int n = static_cast<int>(arr.size());
  std::vector<int> buffer(n);

  for (int byte_idx = 0; byte_idx < 4; ++byte_idx) {
    std::array<int, 256> count{};
    bool is_last_byte = (byte_idx == 3);

    for (int x : arr) {
      unsigned char b = (static_cast<unsigned int>(x) >> (byte_idx * 8)) & 0xFF;
      if (is_last_byte) {
        b ^= 0x80;
      }
      count.at(b)++;
    }

    std::array<size_t, 256> offsets{};
    offsets.at(0) = 0;
    for (int i = 1; i < 256; ++i) {
      offsets.at(i) = offsets.at(i - 1) + count.at(i - 1);
    }

    for (int x : arr) {
      unsigned char b = (static_cast<unsigned int>(x) >> (byte_idx * 8)) & 0xFF;
      if (is_last_byte) {
        b ^= 0x80;
      }
      buffer.at(offsets.at(b)++) = x;
    }
    arr = buffer;
  }
}

void LevonychevIRadixBatcherSortALL::NetworkMergeAndSplit(std::vector<int> &local_data, int partner, bool keep_low) {
  int local_size = static_cast<int>(local_data.size());
  int partner_size = 0;

  MPI_Sendrecv(&local_size, 1, MPI_INT, partner, 0, &partner_size, 1, MPI_INT, partner, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

  std::vector<int> partner_data(partner_size);
  MPI_Sendrecv(local_data.data(), local_size, MPI_INT, partner, 1, partner_data.data(), partner_size, MPI_INT, partner,
               1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  std::vector<int> merged;
  merged.reserve(local_size + partner_size);
  std::ranges::merge(local_data, partner_data, std::back_inserter(merged));

  if (keep_low) {
    local_data.assign(merged.begin(), merged.begin() + local_size);
  } else {
    local_data.assign(merged.end() - local_size, merged.end());
  }
}

void LevonychevIRadixBatcherSortALL::CalculateDistribution(int total_n, int size, std::vector<int> &counts,
                                                           std::vector<int> &displs) {
  int base_size = total_n / size;
  int extra = total_n % size;
  for (int i = 0; i < size; ++i) {
    counts[i] = base_size + (i < extra ? 1 : 0);
    displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
  }
}

void LevonychevIRadixBatcherSortALL::LocalSortPhase(std::vector<int> &local_data) {
  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  int loc_n = static_cast<int>(local_data.size());

  std::vector<std::vector<int>> blocks(num_threads);
  int current_pos = 0;
  for (int i = 0; i < num_threads; ++i) {
    int b_size = (loc_n / num_threads) + (i < (loc_n % num_threads) ? 1 : 0);
    if (b_size > 0) {
      blocks[i].assign(local_data.begin() + current_pos, local_data.begin() + current_pos + b_size);
      current_pos += b_size;
    }
  }

  std::vector<std::future<void>> sort_futures;
  for (auto &b : blocks) {
    if (!b.empty()) {
      sort_futures.push_back(std::async(std::launch::async, [&b]() { LocalRadixSort(b); }));
    }
  }
  for (auto &f : sort_futures) {
    f.wait();
  }

  LocalBatcherMerge(blocks);

  local_data.clear();
  for (const auto &b : blocks) {
    local_data.insert(local_data.end(), b.begin(), b.end());
  }
}

void LevonychevIRadixBatcherSortALL::CompareAndMergeBlocks(std::vector<int> &b1, std::vector<int> &b2) {
  if (b1.empty() || b2.empty()) {
    return;
  }

  std::vector<int> merged;
  merged.reserve(b1.size() + b2.size());
  std::ranges::merge(b1, b2, std::back_inserter(merged));

  auto mid = b1.size();
  auto mid_diff = static_cast<std::ptrdiff_t>(mid);

  b1.assign(merged.begin(), merged.begin() + mid_diff);
  b2.assign(merged.begin() + mid_diff, merged.end());
}

void LevonychevIRadixBatcherSortALL::BatcherStep(std::vector<std::vector<int>> &blocks, int pr, int k) {
  int n_blocks = static_cast<int>(blocks.size());
  std::vector<std::future<void>> m_futures;

  for (int j = k % pr; j <= n_blocks - 1 - k; j += 2 * k) {
    for (int i = 0; i < std::min(k, n_blocks - j - k); ++i) {
      int i1 = j + i;
      int i2 = j + i + k;

      if ((i1 / (pr * 2)) == (i2 / (pr * 2))) {
        m_futures.push_back(
            std::async(std::launch::async, [&blocks, i1, i2]() { CompareAndMergeBlocks(blocks[i1], blocks[i2]); }));
      }
    }
  }

  for (auto &f : m_futures) {
    f.wait();
  }
}

void LevonychevIRadixBatcherSortALL::LocalBatcherMerge(std::vector<std::vector<int>> &blocks) {
  int n_blocks = static_cast<int>(blocks.size());

  for (int pr = 1; pr < n_blocks; pr <<= 1) {
    for (int k = pr; k > 0; k >>= 1) {
      BatcherStep(blocks, pr, k);
    }
  }
}

void LevonychevIRadixBatcherSortALL::GlobalCompareExchange(std::vector<int> &local_data, int rank, int i1, int i2) {
  if (rank == i1) {
    NetworkMergeAndSplit(local_data, i2, true);
  } else if (rank == i2) {
    NetworkMergeAndSplit(local_data, i1, false);
  }
}
void LevonychevIRadixBatcherSortALL::GlobalBatcherStep(std::vector<int> &local_data, int rank, int size, int p, int k) {
  for (int j = k % p; j <= size - 1 - k; j += 2 * k) {
    for (int i = 0; i < std::min(k, size - j - k); ++i) {
      int idx1 = j + i;
      int idx2 = j + i + k;

      if ((idx1 / (p * 2)) == (idx2 / (p * 2))) {
        GlobalCompareExchange(local_data, rank, idx1, idx2);
      }
    }
  }
}

void LevonychevIRadixBatcherSortALL::GlobalSortPhase(std::vector<int> &local_data, int rank, int size) {
  for (int pr = 1; pr < size; pr <<= 1) {
    for (int k = pr; k > 0; k >>= 1) {
      GlobalBatcherStep(local_data, rank, size, pr, k);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

bool LevonychevIRadixBatcherSortALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const std::vector<int> &input = GetInput();
  int total_n = static_cast<int>(input.size());

  std::vector<int> send_counts(size);
  std::vector<int> displs(size);
  CalculateDistribution(total_n, size, send_counts, displs);

  std::vector<int> local_data(send_counts[rank]);
  MPI_Scatterv(input.data(), send_counts.data(), displs.data(), MPI_INT, local_data.data(), send_counts[rank], MPI_INT,
               0, MPI_COMM_WORLD);
  LocalSortPhase(local_data);

  GlobalSortPhase(local_data, rank, size);

  std::vector<int> result;
  if (rank == 0) {
    result.resize(total_n);
  }

  MPI_Gatherv(local_data.data(), static_cast<int>(local_data.size()), MPI_INT, result.data(), send_counts.data(),
              displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    GetOutput() = result;
  }

  return true;
}

bool LevonychevIRadixBatcherSortALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return !GetInput().empty();
  }
  return true;
}
bool LevonychevIRadixBatcherSortALL::PreProcessingImpl() {
  return true;
}
bool LevonychevIRadixBatcherSortALL::PostProcessingImpl() {
  return true;
}

}  // namespace levonychev_i_radix_batcher_sort
