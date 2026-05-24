#include "chetverikova_e_shell_sort_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

ChetverikovaEShellSortSimpleMergeALL::ChetverikovaEShellSortSimpleMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    GetInput() = in;
  } else {
    GetInput().clear();
  }

  GetOutput().clear();
}

bool ChetverikovaEShellSortSimpleMergeALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    return !GetInput().empty();
  }

  return true;
}

bool ChetverikovaEShellSortSimpleMergeALL::PreProcessingImpl() {
  return true;
}

void ChetverikovaEShellSortSimpleMergeALL::ShellSort(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  size_t n = data.size();

  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; ++i) {
      int temp = data[i];
      size_t j = i;

      while (j >= gap && data[j - gap] > temp) {
        data[j] = data[j - gap];
        j -= gap;
      }

      data[j] = temp;
    }
  }
}

std::vector<int> ChetverikovaEShellSortSimpleMergeALL::MergeTwoSortedVectors(const std::vector<int> &a,
                                                                             const std::vector<int> &b) {
  std::vector<int> result(a.size() + b.size());

  std::ranges::merge(a, b, result.begin());

  return result;
}

void ChetverikovaEShellSortSimpleMergeALL::CalculateCountsAndDisplacements(int global_size, int processes_count,
                                                                           std::vector<int> &counts,
                                                                           std::vector<int> &displacements) {
  const int base = global_size / processes_count;
  const int remainder = global_size % processes_count;

  int offset = 0;

  for (int i = 0; i < processes_count; ++i) {
    counts[i] = base + (i < remainder ? 1 : 0);
    displacements[i] = offset;
    offset += counts[i];
  }
}

std::vector<int> ChetverikovaEShellSortSimpleMergeALL::MergeLocalBuffers(std::vector<std::vector<int>> &local_buffers) {
  if (local_buffers.empty()) {
    return {};
  }

  std::vector<int> result = std::move(local_buffers[0]);

  for (size_t i = 1; i < local_buffers.size(); ++i) {
    result = MergeTwoSortedVectors(result, local_buffers[i]);
  }

  return result;
}

bool ChetverikovaEShellSortSimpleMergeALL::RunImpl() {
  int rank = 0;
  int size = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int global_size = 0;

  if (rank == 0) {
    global_size = static_cast<int>(GetInput().size());
  }

  MPI_Bcast(&global_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> sendcounts(size);
  std::vector<int> displs(size);

  if (rank == 0) {
    CalculateCountsAndDisplacements(global_size, size, sendcounts, displs);
  }

  MPI_Bcast(sendcounts.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(displs.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

  int local_size = 0;

  MPI_Scatter(sendcounts.data(), 1, MPI_INT, &local_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> local_data(local_size);

  const int *input_ptr = nullptr;
  if (rank == 0) {
    input_ptr = GetInput().data();
  }

  MPI_Scatterv(input_ptr, sendcounts.data(), displs.data(), MPI_INT, local_data.data(), local_size, MPI_INT, 0,
               MPI_COMM_WORLD);

  const size_t threads = std::max(1, omp_get_max_threads());
  const size_t parts = std::max<size_t>(1, std::min<size_t>(threads, local_data.size()));

  std::vector<std::vector<int>> buffers(parts);
  std::vector<size_t> borders(parts + 1, 0);

  const size_t block = local_data.size() / parts;
  const size_t rem = local_data.size() % parts;

  for (size_t i = 0; i < parts; ++i) {
    borders[i + 1] = borders[i] + block + (i < rem ? 1 : 0);
  }

#pragma omp parallel for default(none) shared(local_data, buffers, borders, parts)
  for (size_t i = 0; i < parts; ++i) {
    using DiffT = std::vector<int>::difference_type;

    auto begin = local_data.begin() + static_cast<DiffT>(borders[i]);
    auto end = local_data.begin() + static_cast<DiffT>(borders[i + 1]);

    std::vector<int> temp(begin, end);

    ShellSort(temp);
    buffers[i] = std::move(temp);
  }

  std::vector<int> local_sorted = MergeLocalBuffers(buffers);

  int local_size_sorted = static_cast<int>(local_sorted.size());

  std::vector<int> recvcounts(size);
  std::vector<int> recvdispls(size);

  MPI_Gather(&local_size_sorted, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  int total = 0;

  if (rank == 0) {
    for (int i = 0; i < size; ++i) {
      recvdispls[i] = total;
      total += recvcounts[i];
    }
    GetOutput().resize(total);
  }

  MPI_Gatherv(local_sorted.data(), local_size_sorted, MPI_INT, (rank == 0 ? GetOutput().data() : nullptr),
              recvcounts.data(), recvdispls.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::vector<std::vector<int>> chunks(size);

    for (int i = 0; i < size; ++i) {
      chunks[i] =
          std::vector<int>(GetOutput().begin() + recvdispls[i], GetOutput().begin() + recvdispls[i] + recvcounts[i]);
    }

    GetOutput() = MergeLocalBuffers(chunks);
  }

  return true;
}

bool ChetverikovaEShellSortSimpleMergeALL::PostProcessingImpl() {
  int rank = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int output_size = 0;

  if (rank == 0) {
    output_size = static_cast<int>(GetOutput().size());
  }

  MPI_Bcast(&output_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    GetOutput().resize(output_size);
  }

  MPI_Bcast(GetOutput().data(), output_size, MPI_INT, 0, MPI_COMM_WORLD);

  return true;
}

}  // namespace chetverikova_e_shell_sort_simple_merge
