#include "zaharov_g_linear_contrast_stretch/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"
#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"

namespace zaharov_g_linear_contrast_stretch {

namespace {

struct MinMax {
  int min;
  int max;
};

std::size_t GetThreadCount(std::size_t size) {
  const auto requested_threads = static_cast<std::size_t>(std::max(1, ppc::util::GetNumThreads()));
  return std::max<std::size_t>(1, std::min(size, requested_threads));
}

std::pair<std::size_t, std::size_t> GetRankRange(std::size_t size) {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto rank_id = static_cast<std::size_t>(rank);
  const auto ranks_count = static_cast<std::size_t>(world_size);
  return {size * rank_id / ranks_count, size * (rank_id + 1) / ranks_count};
}

MinMax FindLocalMinMaxStl(const InType &input, std::size_t begin, std::size_t end) {
  if (begin == end) {
    return {.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};
  }

  const std::size_t range_size = end - begin;
  const std::size_t thread_count = GetThreadCount(range_size);
  std::vector<MinMax> local_minmax(
      thread_count, {.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()});
  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  for (std::size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    const std::size_t thread_begin = begin + (range_size * thread_index / thread_count);
    const std::size_t thread_end = begin + (range_size * (thread_index + 1) / thread_count);
    threads.emplace_back([&input, &local_minmax, thread_begin, thread_end, thread_index]() {
      MinMax current{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};
      for (std::size_t i = thread_begin; i < thread_end; ++i) {
        const int value = static_cast<int>(input[i]);
        current.min = std::min(current.min, value);
        current.max = std::max(current.max, value);
      }
      local_minmax[thread_index] = current;
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  MinMax result{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};
  for (const auto &current : local_minmax) {
    result.min = std::min(result.min, current.min);
    result.max = std::max(result.max, current.max);
  }
  return result;
}

MinMax FindGlobalMinMax(const InType &input) {
  const auto [begin, end] = GetRankRange(input.size());
  const MinMax local = FindLocalMinMaxStl(input, begin, end);

  MinMax global{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};
  MPI_Allreduce(&local.min, &global.min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local.max, &global.max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return global;
}

void StretchOmpRange(const InType &input, OutType &output, std::size_t begin, std::size_t end, int min_el, int denom) {
  const auto signed_begin = static_cast<std::int64_t>(begin);
  const auto signed_end = static_cast<std::int64_t>(end);
#pragma omp parallel for default(none) shared(input, output, min_el, denom, signed_begin, signed_end) \
    num_threads(ppc::util::GetNumThreads())
  for (std::int64_t i = signed_begin; i < signed_end; ++i) {
    const int value =
        (static_cast<int>(input[static_cast<std::size_t>(i)]) - min_el) * std::numeric_limits<uint8_t>::max() / denom;
    output[static_cast<std::size_t>(i)] = static_cast<uint8_t>(std::clamp(value, 0, 255));
  }
}

void StretchStlRange(const InType &input, OutType &output, std::size_t begin, std::size_t end, int min_el, int denom) {
  if (begin == end) {
    return;
  }

  const std::size_t range_size = end - begin;
  const std::size_t thread_count = GetThreadCount(range_size);
  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  for (std::size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    const std::size_t thread_begin = begin + (range_size * thread_index / thread_count);
    const std::size_t thread_end = begin + (range_size * (thread_index + 1) / thread_count);
    threads.emplace_back([&input, &output, thread_begin, thread_end, min_el, denom]() {
      for (std::size_t i = thread_begin; i < thread_end; ++i) {
        const int value = (static_cast<int>(input[i]) - min_el) * std::numeric_limits<uint8_t>::max() / denom;
        output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

void StretchTbbRange(const InType &input, OutType &output, std::size_t begin, std::size_t end, int min_el, int denom) {
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(begin, end),
                            [&input, &output, min_el, denom](const oneapi::tbb::blocked_range<std::size_t> &range) {
    for (std::size_t i = range.begin(); i != range.end(); ++i) {
      const int value = (static_cast<int>(input[i]) - min_el) * std::numeric_limits<uint8_t>::max() / denom;
      output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
    }
  });
}

void CopyTbb(const InType &input, OutType &output) {
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, input.size()),
                            [&input, &output](const oneapi::tbb::blocked_range<std::size_t> &range) {
    std::copy(input.begin() + static_cast<std::ptrdiff_t>(range.begin()),
              input.begin() + static_cast<std::ptrdiff_t>(range.end()),
              output.begin() + static_cast<std::ptrdiff_t>(range.begin()));
  });
}

void StretchImage(const InType &input, OutType &output, const MinMax &minmax) {
  if (minmax.max <= minmax.min) {
    CopyTbb(input, output);
    return;
  }

  const int denom = minmax.max - minmax.min;
  const std::size_t first_border = input.size() / 3;
  const std::size_t second_border = (2 * input.size()) / 3;

  StretchOmpRange(input, output, 0, first_border, minmax.min, denom);
  StretchStlRange(input, output, first_border, second_border, minmax.min, denom);
  StretchTbbRange(input, output, second_border, input.size(), minmax.min, denom);
}

}  // namespace

ZaharovGLinContrStrALL::ZaharovGLinContrStrALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZaharovGLinContrStrALL::ValidationImpl() {
  return !GetInput().empty();
}

bool ZaharovGLinContrStrALL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool ZaharovGLinContrStrALL::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  const MinMax minmax = FindGlobalMinMax(input);
  StretchImage(input, output, minmax);

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool ZaharovGLinContrStrALL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace zaharov_g_linear_contrast_stretch
