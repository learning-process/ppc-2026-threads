#include "batkov_f_contrast_enh_lin_hist_stretch/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "batkov_f_contrast_enh_lin_hist_stretch/common/include/common.hpp"
#include "util/include/util.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch {

namespace {

void BuildScatterLayout(size_t comm_size, size_t image_size, std::vector<int> &sendcounts, std::vector<int> &displs) {
  sendcounts.resize(comm_size);
  displs.resize(comm_size);

  size_t offset = 0;
  const size_t base = image_size / comm_size;
  const size_t rem = image_size % comm_size;

  for (size_t rank = 0; rank < comm_size; ++rank) {
    sendcounts[rank] = static_cast<int>(base + (rank < rem ? 1 : 0));
    displs[rank] = static_cast<int>(offset);

    offset += sendcounts[rank];
  }
}

std::pair<uint8_t, uint8_t> FindMinMaxParallel(const InType &chunk, size_t num_threads) {
  const size_t chunk_size = chunk.size();
  const size_t block = chunk_size / num_threads;

  std::vector<uint8_t> mins(num_threads, std::numeric_limits<uint8_t>::max());
  std::vector<uint8_t> maxs(num_threads, std::numeric_limits<uint8_t>::min());

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (size_t thread_index = 0; thread_index < num_threads; ++thread_index) {
    const size_t begin = thread_index * block;
    const size_t end = (thread_index == num_threads - 1) ? chunk_size : begin + block;

    threads.emplace_back([&, thread_index, begin, end]() {
      uint8_t local_min = std::numeric_limits<uint8_t>::max();
      uint8_t local_max = std::numeric_limits<uint8_t>::min();

      for (size_t i = begin; i < end; ++i) {
        local_min = std::min(local_min, chunk[i]);
        local_max = std::max(local_max, chunk[i]);
      }

      mins[thread_index] = local_min;
      maxs[thread_index] = local_max;
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  return {*std::ranges::min_element(mins), *std::ranges::max_element(maxs)};
}

std::pair<uint8_t, uint8_t> FindLocalMinMax(const InType &chunk, size_t parallel_threshold, size_t num_threads) {
  if (chunk.empty()) {
    return {std::numeric_limits<uint8_t>::max(), std::numeric_limits<uint8_t>::min()};
  }

  if (chunk.size() < parallel_threshold || num_threads <= 1) {
    const auto [min_it, max_it] = std::ranges::minmax_element(chunk);
    return {*min_it, *max_it};
  }

  return FindMinMaxParallel(chunk, num_threads);
}

std::array<uint8_t, 256> BuildStretchLut(float coeff_a, float coeff_b) {
  std::array<uint8_t, 256> lut{};
  for (size_t pixel = 0; pixel < 256; ++pixel) {
    lut.at(pixel) = static_cast<uint8_t>(std::clamp((coeff_a * static_cast<float>(pixel)) + coeff_b, 0.0F, 255.0F));
  }

  return lut;
}

void ApplyStretchParallel(const InType &chunk_in, OutType &chunk_out, size_t num_threads,
                          const std::array<uint8_t, 256> &lut) {
  const size_t n = chunk_in.size();
  chunk_out.resize(n);

  if (n == 0) {
    return;
  }

  if (num_threads <= 1) {
    for (size_t i = 0; i < n; ++i) {
      chunk_out[i] = lut.at(chunk_in[i]);
    }

    return;
  }

  const size_t block_size = n / num_threads;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (size_t thread_index = 0; thread_index < num_threads; ++thread_index) {
    const size_t begin = thread_index * block_size;
    const size_t end = (thread_index == num_threads - 1) ? n : begin + block_size;

    threads.emplace_back([&, begin, end]() {
      for (size_t i = begin; i < end; ++i) {
        chunk_out[i] = lut.at(chunk_in[i]);
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }
}

void CopyChunkParallel(const InType &chunk_in, OutType &chunk_out, size_t num_threads) {
  const size_t n = chunk_in.size();
  chunk_out.resize(n);
  if (n == 0) {
    return;
  }

  if (num_threads <= 1) {
    std::ranges::copy(chunk_in, chunk_out.begin());
    return;
  }

  const size_t block_size = n / num_threads;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (size_t thread_index = 0; thread_index < num_threads; ++thread_index) {
    const size_t begin = thread_index * block_size;
    const size_t end = (thread_index == num_threads - 1) ? n : begin + block_size;

    threads.emplace_back([&, begin, end]() {
      for (size_t i = begin; i < end; ++i) {
        chunk_out[i] = chunk_in[i];
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }
}

}  // namespace

BatkovFContrastEnhLinHistStretchALL::BatkovFContrastEnhLinHistStretchALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size_);

  if (rank_ == 0) {
    GetInput() = in;
  }
}

bool BatkovFContrastEnhLinHistStretchALL::ValidationImpl() {
  if (rank_ == 0) {
    return !GetInput().empty();
  }
  return GetInput().empty();
}

bool BatkovFContrastEnhLinHistStretchALL::PreProcessingImpl() {
  uint64_t input_size_u64 = 0;
  if (rank_ == 0) {
    input_size_u64 = static_cast<uint64_t>(GetInput().size());
  }
  MPI_Bcast(&input_size_u64, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (input_size_u64 > static_cast<uint64_t>(std::numeric_limits<std::ptrdiff_t>::max())) {
    return false;
  }

  image_size_ = static_cast<size_t>(input_size_u64);
  GetOutput().resize(image_size_);
  return image_size_ > 0;
}

bool BatkovFContrastEnhLinHistStretchALL::RunImpl() {
  const auto rank = static_cast<size_t>(rank_);
  const auto comm_size = static_cast<size_t>(comm_size_);

  std::vector<int> sendcounts;
  std::vector<int> displs;
  BuildScatterLayout(comm_size, image_size_, sendcounts, displs);

  const int recvcount = sendcounts[rank];
  std::vector<uint8_t> local_input(static_cast<size_t>(recvcount));
  std::vector<uint8_t> local_output;

  const uint8_t *sendbuf = (rank == 0) ? GetInput().data() : nullptr;

  MPI_Scatterv(sendbuf, sendcounts.data(), displs.data(), MPI_UNSIGNED_CHAR, local_input.data(), recvcount,
               MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  constexpr size_t kParallelMinMaxThreshold = 100000;
  const size_t num_threads = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));

  auto [local_min, local_max] = FindLocalMinMax(local_input, kParallelMinMaxThreshold, num_threads);

  unsigned char global_min = 0;
  unsigned char global_max = 0;
  MPI_Allreduce(&local_min, &global_min, 1, MPI_UNSIGNED_CHAR, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_max, &global_max, 1, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD);

  if (global_min == global_max) {
    CopyChunkParallel(local_input, local_output, num_threads);
  } else {
    const float coeff_a = 255.0F / static_cast<float>(global_max - global_min);
    const float coeff_b = -coeff_a * static_cast<float>(global_min);
    const auto lut = BuildStretchLut(coeff_a, coeff_b);
    ApplyStretchParallel(local_input, local_output, num_threads, lut);
  }

  auto &output = GetOutput();
  uint8_t *gather_recv = (rank == 0) ? output.data() : nullptr;
  MPI_Gatherv(local_output.data(), recvcount, MPI_UNSIGNED_CHAR, gather_recv, sendcounts.data(), displs.data(),
              MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  MPI_Bcast(output.data(), static_cast<int>(image_size_), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  return true;
}

bool BatkovFContrastEnhLinHistStretchALL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace batkov_f_contrast_enh_lin_hist_stretch
