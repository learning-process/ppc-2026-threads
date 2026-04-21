#include "mityaeva_radix/all/include/ops_all.hpp"

#include <mpi.h>

#include <cstddef>
#include <cstring>
#include <utility>
#include <vector>

#include "mityaeva_radix/common/include/common.hpp"
#include "mityaeva_radix/omp/include/sorter_omp.hpp"

namespace mityaeva_radix {

namespace {

void ComputeChunkParams(size_t total_size, int mpi_size, std::vector<size_t> &chunk_sizes,
                        std::vector<size_t> &offsets) {
  size_t base_chunk = total_size / mpi_size;
  size_t remainder = total_size % mpi_size;

  for (int i = 0; i < mpi_size; ++i) {
    chunk_sizes[i] = base_chunk + (std::cmp_less(i, remainder) ? 1 : 0);
    offsets[i] = (i == 0) ? 0 : offsets[i - 1] + chunk_sizes[i - 1];
  }
}

void ScatterData(const std::vector<double> &array, std::vector<double> &local_data,
                 const std::vector<size_t> &chunk_sizes, const std::vector<size_t> &offsets) {
  int mpi_size = static_cast<int>(chunk_sizes.size());
  std::vector<int> send_counts(mpi_size);
  std::vector<int> send_displs(mpi_size);

  for (int i = 0; i < mpi_size; ++i) {
    send_counts[i] = static_cast<int>(chunk_sizes[i]);
    send_displs[i] = static_cast<int>(offsets[i]);
  }

  MPI_Scatterv(array.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE, local_data.data(),
               static_cast<int>(local_data.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

std::vector<double> MergeTwoSorted(const std::vector<double> &left, const std::vector<double> &right) {
  std::vector<double> result(left.size() + right.size());
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

void ExchangeAndMerge(int partner, std::vector<double> &merged_data) {
  size_t my_size = merged_data.size();
  size_t partner_size = 0;

  MPI_Sendrecv(&my_size, 1, MPI_UNSIGNED_LONG, partner, 0, &partner_size, 1, MPI_UNSIGNED_LONG, partner, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  std::vector<double> partner_data(partner_size);
  MPI_Sendrecv(merged_data.data(), static_cast<int>(my_size), MPI_DOUBLE, partner, 1, partner_data.data(),
               static_cast<int>(partner_size), MPI_DOUBLE, partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  merged_data = MergeTwoSorted(merged_data, partner_data);
}

void ParallelHypercubeMerge(std::vector<double> &merged_data, int mpi_rank, int mpi_size) {
  int step = 1;

  while (step < mpi_size) {
    int partner = mpi_rank ^ step;

    if (partner < mpi_size) {
      ExchangeAndMerge(partner, merged_data);
    }

    step <<= 1;
  }
}

}  // namespace

MityaevaRadixAll::MityaevaRadixAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MityaevaRadixAll::ValidationImpl() {
  return !GetInput().empty();
}

bool MityaevaRadixAll::PreProcessingImpl() {
  return true;
}

bool MityaevaRadixAll::RunImpl() {
  int mpi_rank = 0;
  int mpi_size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  auto &array = GetInput();
  size_t total_size = array.size();

  if (total_size == 0) {
    GetOutput() = array;
    return true;
  }

  std::vector<size_t> chunk_sizes(mpi_size);
  std::vector<size_t> offsets(mpi_size);
  ComputeChunkParams(total_size, mpi_size, chunk_sizes, offsets);

  std::vector<double> local_data(chunk_sizes[mpi_rank]);
  ScatterData(array, local_data, chunk_sizes, offsets);

  SorterOmp::Sort(local_data);

  if (mpi_size == 1) {
    GetOutput() = local_data;
    return true;
  }

  std::vector<double> merged_data = std::move(local_data);
  ParallelHypercubeMerge(merged_data, mpi_rank, mpi_size);

  if (mpi_rank == 0) {
    GetOutput() = merged_data;
  } else {
    GetOutput() = {};
  }

  return true;
}

bool MityaevaRadixAll::PostProcessingImpl() {
  return true;
}

}  // namespace mityaeva_radix
