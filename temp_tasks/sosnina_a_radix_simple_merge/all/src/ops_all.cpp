#include "sosnina_a_radix_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "sosnina_a_radix_simple_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace sosnina_a_radix_simple_merge {

namespace {

constexpr int kRadixBits = 8;
constexpr int kRadixSize = 1 << kRadixBits;
constexpr int kNumPasses = sizeof(int) / sizeof(uint8_t);
constexpr uint32_t kSignFlip = 0x80000000U;
constexpr size_t kMinElementsPerPart = 4096;
constexpr size_t kMinElementsPerPartSmall = 32768;
constexpr size_t kSmallArrayThreshold = 1'000'000;
constexpr size_t kLargeArrayThreshold = 20'000'000;

void RadixSortLSD(std::vector<int> &data, std::vector<int> &buffer) {
  size_t idx = 0;
  for (int elem : data) {
    buffer[idx++] = static_cast<int>(static_cast<uint32_t>(elem) ^ kSignFlip);
  }
  std::swap(data, buffer);

  for (int pass = 0; pass < kNumPasses; ++pass) {
    std::array<int, kRadixSize + 1> count{};

    for (auto elem : data) {
      auto digit = static_cast<uint8_t>((static_cast<uint32_t>(elem) >> (pass * kRadixBits)) & 0xFF);
      ++count.at(static_cast<size_t>(digit) + 1U);
    }

    for (int i = 1; i <= kRadixSize; ++i) {
      const auto ui = static_cast<size_t>(i);
      count.at(ui) += count.at(ui - 1U);
    }

    for (auto elem : data) {
      auto digit = static_cast<uint8_t>((static_cast<uint32_t>(elem) >> (pass * kRadixBits)) & 0xFF);
      const auto di = static_cast<size_t>(digit);
      const int write_pos = count.at(di)++;
      buffer[static_cast<size_t>(write_pos)] = elem;
    }

    std::swap(data, buffer);
  }

  for (int &elem : data) {
    elem = static_cast<int>(static_cast<uint32_t>(elem) ^ kSignFlip);
  }
}

void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result) {
  std::ranges::merge(left, right, result.begin());
}

template <typename F>
void ParallelForRange(size_t begin, size_t end, int num_threads, F &&fn) {
  if (begin >= end) {
    return;
  }
  std::decay_t<F> func{std::forward<F>(fn)};
  num_threads = std::max(1, std::min(num_threads, static_cast<int>(end - begin)));
  const size_t n = end - begin;
  const size_t chunk = (n + static_cast<size_t>(num_threads) - 1) / static_cast<size_t>(num_threads);
  std::vector<std::thread> threads;
  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    const size_t lo = begin + (static_cast<size_t>(thread_idx) * chunk);
    if (lo >= end) {
      break;
    }
    const size_t hi = std::min(end, lo + chunk);
    threads.emplace_back([lo, hi, &func]() {
      for (size_t i = lo; i < hi; ++i) {
        func(i);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
}

/// Локальная параллельная radix + merge (как STL): потоки внутри MPI-процесса.
void SortLocalStlParallel(std::vector<int> &data) {
  if (data.size() <= 1) {
    return;
  }

  const int num_threads = ppc::util::GetNumThreads();
  const size_t per_thread_floor = data.size() >= kLargeArrayThreshold
                                      ? (data.size() / static_cast<size_t>(std::max(1, num_threads)))
                                      : (data.size() / static_cast<size_t>(std::max(1, 2 * num_threads)));
  size_t min_chunk_base = kMinElementsPerPart;
  if (data.size() < kMinElementsPerPartSmall) {
    min_chunk_base = std::max(size_t{1}, data.size() / static_cast<size_t>(std::max(1, 2 * num_threads)));
  } else if (data.size() < kSmallArrayThreshold) {
    min_chunk_base = kMinElementsPerPartSmall;
  }
  const size_t min_chunk = std::max(min_chunk_base, per_thread_floor);
  const int max_parts_by_grain = std::max(1, static_cast<int>(data.size() / min_chunk));
  const int num_parts = std::min({num_threads, static_cast<int>(data.size()), max_parts_by_grain});

  if (num_parts <= 1) {
    std::vector<int> buffer(data.size());
    RadixSortLSD(data, buffer);
    return;
  }

  std::vector<std::vector<int>> parts(static_cast<size_t>(num_parts));
  const size_t base_size = data.size() / static_cast<size_t>(num_parts);
  const size_t remainder = data.size() % static_cast<size_t>(num_parts);
  size_t pos = 0;

  for (int i = 0; i < num_parts; ++i) {
    const size_t part_size = base_size + (std::cmp_less(i, remainder) ? 1U : 0U);
    parts[static_cast<size_t>(i)].assign(data.begin() + static_cast<std::ptrdiff_t>(pos),
                                         data.begin() + static_cast<std::ptrdiff_t>(pos + part_size));
    pos += part_size;
  }

  ParallelForRange(0, static_cast<size_t>(num_parts), num_threads, [&](size_t i) {
    auto &part = parts[i];
    std::vector<int> buffer(part.size());
    RadixSortLSD(part, buffer);
  });

  std::vector<std::vector<int>> current = std::move(parts);
  while (current.size() > 1) {
    const size_t half = (current.size() + 1) / 2;
    std::vector<std::vector<int>> next(half);

    const size_t pair_count = current.size() / 2;
    ParallelForRange(0, pair_count, num_threads, [&](size_t idx) {
      auto &left = current[2 * idx];
      auto &right = current[(2 * idx) + 1];
      next[idx].resize(left.size() + right.size());
      SimpleMerge(left, right, next[idx]);
      std::vector<int>().swap(left);
      std::vector<int>().swap(right);
    });
    if (current.size() % 2 == 1) {
      next[half - 1] = std::move(current.back());
    }
    current = std::move(next);
  }

  data = std::move(current[0]);
}

void ComputeChunkParams(size_t total_size, int mpi_size, std::vector<size_t> &chunk_sizes,
                        std::vector<size_t> &offsets) {
  const size_t base_chunk = total_size / static_cast<size_t>(mpi_size);
  const size_t remainder = total_size % static_cast<size_t>(mpi_size);

  for (int i = 0; i < mpi_size; ++i) {
    chunk_sizes[static_cast<size_t>(i)] = base_chunk + (std::cmp_less(i, remainder) ? 1U : 0U);
    offsets[static_cast<size_t>(i)] =
        (i == 0) ? 0 : offsets[static_cast<size_t>(i - 1)] + chunk_sizes[static_cast<size_t>(i - 1)];
  }
}

void ScatterData(const std::vector<int> &array, std::vector<int> &local_data, const std::vector<size_t> &chunk_sizes,
                 const std::vector<size_t> &offsets) {
  const int mpi_size = static_cast<int>(chunk_sizes.size());
  std::vector<int> send_counts(mpi_size);
  std::vector<int> send_displs(mpi_size);

  for (int i = 0; i < mpi_size; ++i) {
    send_counts[i] = static_cast<int>(chunk_sizes[static_cast<size_t>(i)]);
    send_displs[i] = static_cast<int>(offsets[static_cast<size_t>(i)]);
  }

  MPI_Scatterv(array.data(), send_counts.data(), send_displs.data(), MPI_INT, local_data.data(),
               static_cast<int>(local_data.size()), MPI_INT, 0, MPI_COMM_WORLD);
}

std::vector<int> MergeTwoSorted(const std::vector<int> &left, const std::vector<int> &right) {
  std::vector<int> result(left.size() + right.size());
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result[k++] = left[i++];
    } else {
      result[k++] = right[j++];
    }
  }

  while (i < left.size()) {
    result[k++] = left[i++];
  }

  while (j < right.size()) {
    result[k++] = right[j++];
  }

  return result;
}

void ExchangeAndMerge(int partner, std::vector<int> &merged_data) {
  const size_t my_size = merged_data.size();
  size_t partner_size = 0;

  MPI_Sendrecv(&my_size, 1, MPI_UNSIGNED_LONG, partner, 0, &partner_size, 1, MPI_UNSIGNED_LONG, partner, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  std::vector<int> partner_data(partner_size);
  MPI_Sendrecv(merged_data.data(), static_cast<int>(my_size), MPI_INT, partner, 1, partner_data.data(),
               static_cast<int>(partner_size), MPI_INT, partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  merged_data = MergeTwoSorted(merged_data, partner_data);
}

void ParallelHypercubeMerge(std::vector<int> &merged_data, int mpi_rank, int mpi_size) {
  int step = 1;

  while (step < mpi_size) {
    const int partner = mpi_rank ^ step;

    if (partner < mpi_size) {
      ExchangeAndMerge(partner, merged_data);
    }

    step <<= 1;
  }
}

void BcastSortedVector(std::vector<int> &data, int mpi_rank) {
  size_t n = data.size();
  MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  if (mpi_rank != 0) {
    data.resize(n);
  }
  if (n > 0U) {
    MPI_Bcast(data.data(), static_cast<int>(n), MPI_INT, 0, MPI_COMM_WORLD);
  }
}

}  // namespace

SosninaATestTaskALL::SosninaATestTaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool SosninaATestTaskALL::ValidationImpl() {
  return !GetInput().empty();
}

bool SosninaATestTaskALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SosninaATestTaskALL::RunImpl() {
  int mpi_rank = 0;
  int mpi_size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  std::vector<int> &data = GetOutput();
  const size_t total_size = data.size();

  if (total_size <= 1) {
    return true;
  }

  std::vector<size_t> chunk_sizes(static_cast<size_t>(mpi_size));
  std::vector<size_t> offsets(static_cast<size_t>(mpi_size));
  ComputeChunkParams(total_size, mpi_size, chunk_sizes, offsets);

  std::vector<int> local_data(chunk_sizes[static_cast<size_t>(mpi_rank)]);
  ScatterData(data, local_data, chunk_sizes, offsets);

  SortLocalStlParallel(local_data);

  if (mpi_size == 1) {
    data = std::move(local_data);
    return std::ranges::is_sorted(data);
  }

  std::vector<int> merged_data = std::move(local_data);
  ParallelHypercubeMerge(merged_data, mpi_rank, mpi_size);

  if (mpi_rank == 0) {
    data = std::move(merged_data);
  } else {
    data.clear();
  }

  BcastSortedVector(data, mpi_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  return std::ranges::is_sorted(data);
}

bool SosninaATestTaskALL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace sosnina_a_radix_simple_merge
