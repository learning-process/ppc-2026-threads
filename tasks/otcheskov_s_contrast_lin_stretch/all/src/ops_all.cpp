#include "otcheskov_s_contrast_lin_stretch/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"
#include "util/include/util.hpp"

namespace otcheskov_s_contrast_lin_stretch {

OtcheskovSContrastLinStretchALL::OtcheskovSContrastLinStretchALL(const InType &in) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput().clear();
  if (rank_ == 0) {
    GetInput() = in;
    GetOutput().resize(GetInput().size());
  } else {
    GetInput().clear();
    GetInput().shrink_to_fit();
  }
}

bool OtcheskovSContrastLinStretchALL::ValidationImpl() {
  if (rank_ == 0) {
    is_valid_ = !GetInput().empty();
  }
  MPI_Bcast(&is_valid_, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  return is_valid_;
}

bool OtcheskovSContrastLinStretchALL::PreProcessingImpl() {
  return true;
}

bool OtcheskovSContrastLinStretchALL::RunImpl() {
  if (!is_valid_) {
    return false;
  }

  const InType &input = GetInput();
  OutType &output = GetOutput();
  int global_size = 0;
  if (rank_ == 0) {
    global_size = static_cast<int>(input.size());
  }
  MPI_Bcast(&global_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> counts(size_);
  std::vector<int> displs(size_);
  int base = global_size / size_;
  int rem = global_size % size_;

  for (int i = 0; i < size_; ++i) {
    counts[i] = base + (i < rem ? 1 : 0);
    displs[i] = (i * base) + std::min(i, rem);
  }
  auto local_size = static_cast<size_t>(counts[rank_]);
  std::vector<uint8_t> local_input(local_size);
  std::vector<uint8_t> local_output(local_size);

  MPI_Scatterv(rank_ == 0 ? input.data() : nullptr, counts.data(), displs.data(), MPI_UINT8_T, local_input.data(),
               static_cast<int>(local_size), MPI_UINT8_T, 0, MPI_COMM_WORLD);

  MinMax local = ComputeMinMax(local_input);

  uint8_t global_min{};
  uint8_t global_max{};

  MPI_Allreduce(&local.min, &global_min, 1, MPI_UINT8_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local.max, &global_max, 1, MPI_UINT8_T, MPI_MAX, MPI_COMM_WORLD);

  if (global_min == global_max) {
    CopyInput(local_input, local_output);
  } else {
    int min_i = static_cast<int>(global_min);
    int range = static_cast<int>(global_max) - min_i;

    LinearStretch(local_input, local_output, min_i, range);
  }

  MPI_Gatherv(local_output.data(), static_cast<int>(local_size), MPI_UINT8_T, rank_ == 0 ? output.data() : nullptr,
              counts.data(), displs.data(), MPI_UINT8_T, 0, MPI_COMM_WORLD);
  return true;
}

bool OtcheskovSContrastLinStretchALL::PostProcessingImpl() {
  return true;
}

OtcheskovSContrastLinStretchALL::MinMax OtcheskovSContrastLinStretchALL::ComputeMinMax(const InType &input) {
  if (input.empty()) {
    return MinMax{};
  }
  const size_t size = input.size();
  const size_t num_threads = std::min<size_t>(static_cast<size_t>(ppc::util::GetNumThreads()), size);

  const size_t block = size / num_threads;
  std::vector<MinMax> local(num_threads, {255, 0});

  {
    std::vector<std::jthread> threads;
    threads.reserve(num_threads);
    for (size_t tid = 0; tid < num_threads; ++tid) {
      size_t begin = tid * block;
      size_t end = (tid == num_threads - 1) ? size : begin + block;

      threads.emplace_back([&, tid, begin, end]() {
        auto [min, max] = std::ranges::minmax_element(input.begin() + static_cast<std::ptrdiff_t>(begin),
                                                      input.begin() + static_cast<std::ptrdiff_t>(end));
        local[tid] = {.min = *min, .max = *max};
      });
    }
  }

  MinMax result{.min = 255, .max = 0};
  for (size_t tid = 0; tid < num_threads; ++tid) {
    result.min = std::min(result.min, local[tid].min);
    result.max = std::max(result.max, local[tid].max);
  }
  return result;
}

void OtcheskovSContrastLinStretchALL::CopyInput(const InType &input, OutType &output) {
  if (input.empty()) {
    return;
  }
  const size_t size = input.size();
  const size_t num_threads = std::min<size_t>(static_cast<size_t>(ppc::util::GetNumThreads()), size);
  const size_t block = size / num_threads;

  std::vector<std::jthread> threads;
  threads.reserve(num_threads);
  for (size_t tid = 0; tid < num_threads; ++tid) {
    size_t begin = tid * block;
    size_t end = (tid == num_threads - 1) ? size : begin + block;

    threads.emplace_back([&, begin, end]() {
      for (size_t i = begin; i < end; ++i) {
        output[i] = input[i];
      }
    });
  }
}

void OtcheskovSContrastLinStretchALL::LinearStretch(const InType &input, OutType &output, int min_i, int range) {
  if (input.empty()) {
    return;
  }
  const size_t size = input.size();
  const size_t num_threads = std::min<size_t>(static_cast<size_t>(ppc::util::GetNumThreads()), size);
  const size_t block = size / num_threads;

  std::vector<std::jthread> threads;
  threads.reserve(num_threads);
  for (size_t tid = 0; tid < num_threads; ++tid) {
    size_t begin = tid * block;
    size_t end = (tid == num_threads - 1) ? size : begin + block;

    threads.emplace_back([&, begin, end, min_i, range]() {
      for (size_t i = begin; i < end; ++i) {
        int pixel = static_cast<int>(input[i]);
        int value = (pixel - min_i) * 255 / (range);
        value = std::clamp(value, 0, 255);
        output[i] = static_cast<uint8_t>(value);
      }
    });
  }
}

}  // namespace otcheskov_s_contrast_lin_stretch
