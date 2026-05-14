#include "frolova_s_radix_sort_double/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/parallel_sort.h>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <vector>

namespace frolova_s_radix_sort_double {

namespace {

void LocalRadixSort(std::vector<double> &chunk) {
  tbb::parallel_sort(chunk.begin(), chunk.end());
}

}  // namespace

FrolovaSRadixSortDoubleALL::FrolovaSRadixSortDoubleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleALL::ValidationImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleALL::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleALL::RunImpl() {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int total_size = 0;
  if (rank == 0) {
    total_size = static_cast<int>(GetInput().size());
  }
  MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total_size == 0) {
    return true;
  }

  std::vector<int> sendcounts(size), displs(size);
  int rem = total_size % size;
  int offset = 0;
  for (int i = 0; i < size; ++i) {
    sendcounts[i] = total_size / size + (i < rem ? 1 : 0);
    displs[i] = offset;
    offset += sendcounts[i];
  }

  std::vector<double> local_data(sendcounts[rank]);
  double dummy;
  double *recvbuf = (sendcounts[rank] > 0) ? local_data.data() : &dummy;

  MPI_Scatterv(rank == 0 ? GetInput().data() : nullptr, sendcounts.data(), displs.data(), MPI_DOUBLE, recvbuf,
               sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  LocalRadixSort(local_data);
  std::vector<double> gathered;
  if (rank == 0) {
    gathered.resize(total_size);
  }

  double *sendbuf = (sendcounts[rank] > 0) ? local_data.data() : &dummy;
  MPI_Gatherv(sendbuf, sendcounts[rank], MPI_DOUBLE, gathered.data(), sendcounts.data(), displs.data(), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  if (rank == 0) {
    std::vector<double> result;
    if (size > 0 && sendcounts[0] > 0) {
      result.assign(gathered.begin(), gathered.begin() + sendcounts[0]);
    }

    for (int i = 1; i < size; ++i) {
      if (sendcounts[i] == 0) {
        continue;
      }
      std::vector<double> temp(result.size() + sendcounts[i]);
      std::merge(result.begin(), result.end(), gathered.begin() + displs[i],
                 gathered.begin() + displs[i] + sendcounts[i], temp.begin());
      result = std::move(temp);
    }
    GetOutput() = std::move(result);
  }

  return true;
}

bool FrolovaSRadixSortDoubleALL::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
