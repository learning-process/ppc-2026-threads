#include "kamaletdinov_r_bitwise_int/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "kamaletdinov_r_bitwise_int/common/include/common.hpp"

namespace kamaletdinov_r_bitwise_int {

namespace {

void CountingSortByDigitALL(std::vector<int> &arr, int exp) {
  const int n = static_cast<int>(arr.size());
  const int thread_count = omp_get_max_threads();
  std::vector<std::array<int, 10>> local_counts(thread_count);

#pragma omp parallel default(none) shared(arr, exp, local_counts, n) num_threads(thread_count)
  {
    const int tid = omp_get_thread_num();
    auto &current = local_counts.at(tid);
    current.fill(0);

#pragma omp for schedule(static)
    for (int i = 0; i < n; i++) {
      const int digit = (arr.at(i) / exp) % 10;
      current.at(digit)++;
    }
  }

  std::array<int, 10> global_count = {};
  for (int ti = 0; ti < thread_count; ti++) {
    for (int di = 0; di < 10; di++) {
      global_count.at(di) += local_counts.at(ti).at(di);
    }
  }

  std::array<int, 10> global_start = {};
  for (int di = 1; di < 10; di++) {
    global_start.at(di) = global_start.at(di - 1) + global_count.at(di - 1);
  }

  std::vector<std::array<int, 10>> thread_offsets(thread_count);
  for (int di = 0; di < 10; di++) {
    int offset = global_start.at(di);
    for (int ti = 0; ti < thread_count; ti++) {
      thread_offsets.at(ti).at(di) = offset;
      offset += local_counts.at(ti).at(di);
    }
  }

  std::vector<int> output(n);

#pragma omp parallel default(none) shared(arr, exp, output, n, thread_count, thread_offsets) num_threads(thread_count)
  {
    const int tid = omp_get_thread_num();
    auto positions = thread_offsets.at(tid);

#pragma omp for schedule(static)
    for (int i = 0; i < n; i++) {
      const int digit = (arr.at(i) / exp) % 10;
      output.at(positions.at(digit)++) = arr.at(i);
    }
  }

  arr.swap(output);
}

void RadixSortPositiveALL(std::vector<int> &arr) {
  if (arr.empty()) {
    return;
  }

  int max_val = *std::ranges::max_element(arr);
  for (int exp = 1; max_val / exp > 0; exp *= 10) {
    CountingSortByDigitALL(arr, exp);
    if (exp > max_val / 10) {
      break;
    }
  }
}

void LocalBitwiseSort(std::vector<int> &arr) {
  if (arr.size() <= 1) {
    return;
  }

  std::vector<int> neg;
  std::vector<int> pos;
  neg.reserve(arr.size());
  pos.reserve(arr.size());

  for (int val : arr) {
    if (val < 0) {
      neg.push_back(-val);
    } else {
      pos.push_back(val);
    }
  }

  RadixSortPositiveALL(neg);
  RadixSortPositiveALL(pos);

  std::size_t idx = 0;
  for (int i = static_cast<int>(neg.size()) - 1; i >= 0; i--) {
    arr.at(idx++) = -neg.at(i);
  }
  for (int v : pos) {
    arr.at(idx++) = v;
  }
}

std::vector<int> MergeSorted(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> result;
  result.reserve(a.size() + b.size());
  std::size_t i = 0;
  std::size_t j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      result.push_back(a[i++]);
    } else {
      result.push_back(b[j++]);
    }
  }
  while (i < a.size()) {
    result.push_back(a[i++]);
  }
  while (j < b.size()) {
    result.push_back(b[j++]);
  }
  return result;
}

}  // namespace

void BitwiseSortALL(std::vector<int> &arr) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int total = static_cast<int>(arr.size());
  MPI_Bcast(&total, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total <= 1) {
    return;
  }

  std::vector<int> send_counts(size);
  std::vector<int> displs(size);
  for (int i = 0; i < size; i++) {
    send_counts[i] = (total / size) + (i < total % size ? 1 : 0);
    displs[i] = (i == 0) ? 0 : displs[i - 1] + send_counts[i - 1];
  }

  std::vector<int> local_data(send_counts[rank]);
  MPI_Scatterv(arr.data(), send_counts.data(), displs.data(), MPI_INT, local_data.data(), send_counts[rank], MPI_INT, 0,
               MPI_COMM_WORLD);

  LocalBitwiseSort(local_data);

  std::vector<int> recv_counts(size);
  MPI_Gather(&send_counts[rank], 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::vector<int> all_data(total);
    std::vector<int> recv_displs(size);
    for (int i = 1; i < size; i++) {
      recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
    }
    MPI_Gatherv(local_data.data(), send_counts[rank], MPI_INT, all_data.data(), recv_counts.data(), recv_displs.data(),
                MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> merged = {all_data.begin(), all_data.begin() + recv_counts[0]};
    for (int i = 1; i < size; i++) {
      std::vector<int> chunk(all_data.begin() + recv_displs[i], all_data.begin() + recv_displs[i] + recv_counts[i]);
      merged = MergeSorted(merged, chunk);
    }
    arr = merged;
  } else {
    MPI_Gatherv(local_data.data(), send_counts[rank], MPI_INT, nullptr, nullptr, nullptr, MPI_INT, 0, MPI_COMM_WORLD);
    arr.resize(total);
  }

  MPI_Bcast(arr.data(), total, MPI_INT, 0, MPI_COMM_WORLD);
}

KamaletdinovRBitwiseIntALL::KamaletdinovRBitwiseIntALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KamaletdinovRBitwiseIntALL::ValidationImpl() {
  return GetInput() >= 0;
}

bool KamaletdinovRBitwiseIntALL::PreProcessingImpl() {
  int n = GetInput();
  data_.resize(n);
  for (int i = 0; i < n; i++) {
    data_[i] = (n / 2) - i;
  }
  return true;
}

bool KamaletdinovRBitwiseIntALL::RunImpl() {
  BitwiseSortALL(data_);
  return true;
}

bool KamaletdinovRBitwiseIntALL::PostProcessingImpl() {
  bool sorted = std::ranges::is_sorted(data_);
  GetOutput() = sorted ? GetInput() : 0;
  return true;
}

}  // namespace kamaletdinov_r_bitwise_int
