#include "batushin_i_incr_contrast_with_lhs/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "batushin_i_incr_contrast_with_lhs/common/include/common.hpp"

namespace batushin_i_incr_contrast_with_lhs {

BatushinIIncrContrastWithLhsALL::BatushinIIncrContrastWithLhsALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  GetInput() = in;
  GetOutput().resize(in.size());
}

bool BatushinIIncrContrastWithLhsALL::ValidationImpl() {
  return !GetInput().empty();
}

bool BatushinIIncrContrastWithLhsALL::PreProcessingImpl() {
  return true;
}

namespace {

unsigned char NormalizePixel(unsigned char pixel, unsigned char min_val, double scale_factor) {
  double normalized = static_cast<double>(pixel - min_val) * scale_factor;

  normalized = std::floor(normalized + 0.5);

  normalized = std::max(normalized, 0.0);
  normalized = std::min(normalized, 255.0);

  return static_cast<unsigned char>(normalized);
}

unsigned int GetNumThreads() {
  unsigned int threads = std::thread::hardware_concurrency();
  return threads == 0 ? 1 : threads;
}

void NormalizeBlock(const std::vector<unsigned char> &source, std::vector<unsigned char> &destination, size_t start,
                    size_t end, unsigned char min_value, double scale) {
  for (size_t i = start; i < end; ++i) {
    destination[i] = NormalizePixel(source[i], min_value, scale);
  }
}

void NormalizeParallel(const std::vector<unsigned char> &source, std::vector<unsigned char> &destination,
                       unsigned char min_value, double scale) {
  const size_t size = source.size();

  unsigned int threads_count = GetNumThreads();

  size_t chunk_size = std::max(static_cast<size_t>(1), size / threads_count);

  size_t blocks = (size + chunk_size - 1) / chunk_size;

  std::vector<std::thread> threads;

  for (size_t block = 0; block < blocks; ++block) {
    size_t start = block * chunk_size;
    size_t end = std::min(start + chunk_size, size);

    threads.emplace_back([&source, &destination, start, end, min_value, scale] {
      NormalizeBlock(source, destination, start, end, min_value, scale);
    });
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }
}

std::pair<unsigned char, unsigned char> FindLocalMinMax(const std::vector<unsigned char> &data) {
  if (data.empty()) {
    return {0, 0};
  }

  unsigned char min_val = data[0];
  unsigned char max_val = data[0];

  for (size_t i = 1; i < data.size(); ++i) {
    min_val = std::min(min_val, data[i]);
    max_val = std::max(max_val, data[i]);
  }

  return {min_val, max_val};
}

}  // namespace

bool BatushinIIncrContrastWithLhsALL::RunImpl() {
  int rank = 0;
  int size = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const std::vector<unsigned char> &input = GetInput();
  std::vector<unsigned char> &output = GetOutput();

  const int global_size = static_cast<int>(input.size());

  if (rank == 0) {
    output.resize(global_size);
  }

  std::vector<int> send_counts(size);
  std::vector<int> displs(size);

  int base = global_size / size;
  int remainder = global_size % size;

  for (int i = 0; i < size; ++i) {
    send_counts[i] = base + (i < remainder ? 1 : 0);

    displs[i] = (i == 0 ? 0 : displs[i - 1] + send_counts[i - 1]);
  }

  std::vector<unsigned char> local_data(send_counts[rank]);

  MPI_Scatterv(input.data(), send_counts.data(), displs.data(), MPI_UNSIGNED_CHAR, local_data.data(), send_counts[rank],
               MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  auto [local_min, local_max] = FindLocalMinMax(local_data);

  unsigned char global_min = 0;
  unsigned char global_max = 0;

  MPI_Allreduce(&local_min, &global_min, 1, MPI_UNSIGNED_CHAR, MPI_MIN, MPI_COMM_WORLD);

  MPI_Allreduce(&local_max, &global_max, 1, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD);

  std::vector<unsigned char> local_result(local_data.size());

  if (global_min == global_max) {
    local_result = local_data;
  } else {
    double scale = 255.0 / static_cast<double>(global_max - global_min);

    NormalizeParallel(local_data, local_result, global_min, scale);
  }

  MPI_Gatherv(local_result.data(), send_counts[rank], MPI_UNSIGNED_CHAR, output.data(), send_counts.data(),
              displs.data(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    output.resize(global_size);
  }

  MPI_Bcast(output.data(), global_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  return true;
}

bool BatushinIIncrContrastWithLhsALL::PostProcessingImpl() {
  return true;
}

}  // namespace batushin_i_incr_contrast_with_lhs
