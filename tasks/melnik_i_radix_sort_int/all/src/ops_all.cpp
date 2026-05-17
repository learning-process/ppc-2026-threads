#include "melnik_i_radix_sort_int/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"
#include "util/include/util.hpp"

namespace melnik_i_radix_sort_int {

namespace {

constexpr int kBitsPerPass = 8;
constexpr std::size_t kBuckets = 1U << kBitsPerPass;

std::pair<std::size_t, std::size_t> ChunkBounds(std::size_t n, int num_ranks, int rank) {
  const auto nz = static_cast<std::size_t>(num_ranks);
  const auto rz = static_cast<std::size_t>(rank);
  const std::size_t base = (n / nz) * rz + std::min(rz, n % nz);
  const std::size_t size = n / nz + (rz < n % nz ? 1U : 0U);
  return {base, base + size};
}

std::pair<std::size_t, std::size_t> ThreadChunkBounds(std::size_t n, int num_threads, int tid) {
  return ChunkBounds(n, num_threads, tid);
}

}  // namespace

MelnikIRadixSortIntALL::MelnikIRadixSortIntALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MelnikIRadixSortIntALL::ValidationImpl() {
  return !GetInput().empty();
}

bool MelnikIRadixSortIntALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return !GetOutput().empty();
}

bool MelnikIRadixSortIntALL::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

void MelnikIRadixSortIntALL::CountingSortByByte(const std::vector<int> &source, std::vector<int> &destination,
                                                std::size_t begin, std::size_t end, std::int64_t exp,
                                                std::int64_t offset) {
  std::array<std::size_t, kBuckets> count{};
  count.fill(0);

  for (std::size_t i = begin; i < end; ++i) {
    const std::int64_t val = static_cast<std::int64_t>(source[i]) + offset;
    const auto bucket = static_cast<std::size_t>((val / exp) % static_cast<std::int64_t>(kBuckets));
    ++count[bucket];
  }

  std::array<std::size_t, kBuckets> pos{};
  pos[0] = begin;
  for (std::size_t b = 1; b < kBuckets; ++b) {
    pos[b] = pos[b - 1U] + count[b - 1U];
  }

  for (std::size_t i = begin; i < end; ++i) {
    const std::int64_t val = static_cast<std::int64_t>(source[i]) + offset;
    const auto bucket = static_cast<std::size_t>((val / exp) % static_cast<std::int64_t>(kBuckets));
    destination[pos[bucket]] = source[i];
    ++pos[bucket];
  }
}

void MelnikIRadixSortIntALL::RadixSortRange(std::vector<int> &data, std::vector<int> &buffer, std::size_t begin,
                                            std::size_t end) {
  if (end - begin <= 1U) {
    return;
  }

  const auto range_begin = data.begin() + static_cast<std::ptrdiff_t>(begin);
  const auto range_end = data.begin() + static_cast<std::ptrdiff_t>(end);
  const auto [min_it, max_it] = std::ranges::minmax_element(range_begin, range_end);
  const auto min_val = static_cast<std::int64_t>(*min_it);
  const auto max_val = static_cast<std::int64_t>(*max_it);
  const std::int64_t offset = (min_val < 0) ? -min_val : 0;
  const std::int64_t max_shifted = max_val + offset;

  std::vector<int> *src = &data;
  std::vector<int> *dst = &buffer;

  for (std::int64_t exp = 1; max_shifted / exp > 0; exp <<= kBitsPerPass) {
    CountingSortByByte(*src, *dst, begin, end, exp, offset);
    std::swap(src, dst);
  }

  // If result ended up in buffer, copy it back to data.
  if (src != &data) {
    std::copy(src->begin() + static_cast<std::ptrdiff_t>(begin), src->begin() + static_cast<std::ptrdiff_t>(end),
              data.begin() + static_cast<std::ptrdiff_t>(begin));
  }
}

void MelnikIRadixSortIntALL::MergeRanges(const std::vector<int> &source, std::vector<int> &destination, Range left,
                                         Range right, std::size_t write_begin) {
  std::size_t li = left.begin;
  std::size_t ri = right.begin;
  std::size_t wi = write_begin;

  while (li < left.end && ri < right.end) {
    if (source[li] <= source[ri]) {
      destination[wi++] = source[li++];
    } else {
      destination[wi++] = source[ri++];
    }
  }

  if (li < left.end) {
    std::copy(source.begin() + static_cast<std::ptrdiff_t>(li), source.begin() + static_cast<std::ptrdiff_t>(left.end),
              destination.begin() + static_cast<std::ptrdiff_t>(wi));
  } else if (ri < right.end) {
    std::copy(source.begin() + static_cast<std::ptrdiff_t>(ri), source.begin() + static_cast<std::ptrdiff_t>(right.end),
              destination.begin() + static_cast<std::ptrdiff_t>(wi));
  }
}

void MelnikIRadixSortIntALL::MergeSortedRangesParallel(std::vector<int> &data, std::vector<int> &buffer,
                                                       const std::vector<Range> &ranges, int num_threads) {
  if (ranges.empty()) {
    return;
  }

  std::vector<Range> cur = ranges;
  std::vector<int> *src = &data;
  std::vector<int> *dst = &buffer;

  while (cur.size() > 1U) {
    const std::size_t pairs = (cur.size() + 1U) / 2U;
    std::vector<Range> next(pairs);

    // Limit concurrency to num_threads, but not more than available pairs.
    const int active = std::min(num_threads, static_cast<int>(pairs));

    if (active <= 1) {
      // Sequential path – avoids thread overhead for tiny workloads.
      for (std::size_t p = 0; p < pairs; ++p) {
        const std::size_t lp = p * 2U;
        const Range left = cur[lp];
        if (lp + 1U >= cur.size()) {
          // Odd range out – copy to dst to keep buffers consistent.
          std::copy(src->begin() + static_cast<std::ptrdiff_t>(left.begin),
                    src->begin() + static_cast<std::ptrdiff_t>(left.end),
                    dst->begin() + static_cast<std::ptrdiff_t>(left.begin));
          next[p] = left;
          continue;
        }
        const Range right = cur[lp + 1U];
        MergeRanges(*src, *dst, left, right, left.begin);
        next[p] = Range{.begin = left.begin, .end = right.end};
      }
    } else {
      // Parallel path – distribute pairs across threads.
      std::vector<std::thread> workers;
      workers.reserve(static_cast<std::size_t>(active));

      for (int tid = 0; tid < active; ++tid) {
        workers.emplace_back([&, tid]() {
          auto [p_begin, p_end] = ThreadChunkBounds(pairs, active, tid);
          for (std::size_t p = p_begin; p < p_end; ++p) {
            const std::size_t lp = p * 2U;
            const Range left = cur[lp];
            if (lp + 1U >= cur.size()) {
              std::copy(src->begin() + static_cast<std::ptrdiff_t>(left.begin),
                        src->begin() + static_cast<std::ptrdiff_t>(left.end),
                        dst->begin() + static_cast<std::ptrdiff_t>(left.begin));
              next[p] = left;
              return;
            }
            const Range right = cur[lp + 1U];
            MergeRanges(*src, *dst, left, right, left.begin);
            next[p] = Range{.begin = left.begin, .end = right.end};
          }
        });
      }
      for (auto &w : workers) {
        w.join();
      }
    }

    cur = std::move(next);
    std::swap(src, dst);
  }

  if (src != &data) {
    data.swap(*src);
  }
}

bool MelnikIRadixSortIntALL::RunImpl() {
  auto &output = GetOutput();
  if (output.empty()) {
    return false;
  }

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  const int num_threads = std::max(1, ppc::util::GetNumThreads());
  const auto total_size = static_cast<int>(output.size());

  std::vector<int> send_counts(static_cast<std::size_t>(num_ranks));
  std::vector<int> displacements(static_cast<std::size_t>(num_ranks));
  if (rank == 0) {
    for (int r = 0; r < num_ranks; ++r) {
      auto [b, e] = ChunkBounds(static_cast<std::size_t>(total_size), num_ranks, r);
      send_counts[static_cast<std::size_t>(r)] = static_cast<int>(e - b);
      displacements[static_cast<std::size_t>(r)] = static_cast<int>(b);
    }
  }

  MPI_Bcast(send_counts.data(), num_ranks, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(displacements.data(), num_ranks, MPI_INT, 0, MPI_COMM_WORLD);

  const int local_count = send_counts[static_cast<std::size_t>(rank)];
  std::vector<int> local(static_cast<std::size_t>(local_count));

  MPI_Scatterv(output.data(), send_counts.data(), displacements.data(), MPI_INT, local.data(), local_count, MPI_INT, 0,
               MPI_COMM_WORLD);

  const int local_threads = num_threads;

  if (local_count > 0) {
    std::vector<int> local_buf(static_cast<std::size_t>(local_count));

    if (local_threads <= 1) {
      RadixSortRange(local, local_buf, 0, static_cast<std::size_t>(local_count));
    } else {
      // Sort each thread's sub-chunk in parallel.
      std::vector<std::thread> workers;
      workers.reserve(static_cast<std::size_t>(local_threads));
      for (int tid = 0; tid < local_threads; ++tid) {
        workers.emplace_back([&, tid]() {
          auto [b, e] = ThreadChunkBounds(static_cast<std::size_t>(local_count), local_threads, tid);
          RadixSortRange(local, local_buf, b, e);
        });
      }
      for (auto &w : workers) {
        w.join();
      }

      // Build ranges for the thread-level merge.
      std::vector<Range> thread_ranges;
      thread_ranges.reserve(static_cast<std::size_t>(local_threads));
      for (int tid = 0; tid < local_threads; ++tid) {
        auto [b, e] = ThreadChunkBounds(static_cast<std::size_t>(local_count), local_threads, tid);
        if (b < e) {
          thread_ranges.push_back(Range{.begin = b, .end = e});
        }
      }

      MergeSortedRangesParallel(local, local_buf, thread_ranges, num_threads);
    }
  }

  MPI_Gatherv(local.data(), local_count, MPI_INT, output.data(), send_counts.data(), displacements.data(), MPI_INT, 0,
              MPI_COMM_WORLD);

  if (rank == 0 && num_ranks > 1) {
    std::vector<int> global_buf(static_cast<std::size_t>(total_size));

    // Build per-rank ranges (already sorted sub-arrays in output).
    std::vector<Range> rank_ranges;
    rank_ranges.reserve(static_cast<std::size_t>(num_ranks));
    for (int r = 0; r < num_ranks; ++r) {
      const auto b = static_cast<std::size_t>(displacements[static_cast<std::size_t>(r)]);
      const auto e = b + static_cast<std::size_t>(send_counts[static_cast<std::size_t>(r)]);
      if (b < e) {
        rank_ranges.push_back(Range{.begin = b, .end = e});
      }
    }

    MergeSortedRangesParallel(output, global_buf, rank_ranges, num_threads);
  }

  MPI_Bcast(output.data(), total_size, MPI_INT, 0, MPI_COMM_WORLD);
  return true;
}

}  // namespace melnik_i_radix_sort_int
