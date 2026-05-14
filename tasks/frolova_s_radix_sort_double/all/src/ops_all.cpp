#include "frolova_s_radix_sort_double/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <vector>

namespace frolova_s_radix_sort_double {

namespace {

void LocalRadixSort(std::vector<double> &chunk) {
  if (chunk.empty()) {
    return;
  }
  const int radix = 256;
  const int num_bits = 8;
  const int num_passes = sizeof(uint64_t);
  std::vector<double> temp(chunk.size());

  for (int pass = 0; pass < num_passes; ++pass) {
    std::vector<int> count(radix, 0);
    for (double val : chunk) {
      auto bits = std::bit_cast<uint64_t>(val);
      count[(bits >> (pass * num_bits)) & 0xFF]++;
    }
    int total = 0;
    for (int i = 0; i < radix; ++i) {
      int old = count[i];
      count[i] = total;
      total += old;
    }
    for (double val : chunk) {
      auto bits = std::bit_cast<uint64_t>(val);
      temp[count[(bits >> (pass * num_bits)) & 0xFF]++] = val;
    }
    chunk.swap(temp);
  }

  std::vector<double> neg, pos;
  for (double val : chunk) {
    if ((std::bit_cast<uint64_t>(val) >> 63) != 0) {
      neg.push_back(val);
    } else {
      pos.push_back(val);
    }
  }
  std::reverse(neg.begin(), neg.end());

  chunk.clear();
  chunk.insert(chunk.end(), neg.begin(), neg.end());
  chunk.insert(chunk.end(), pos.begin(), pos.end());
}

}  // namespace

FrolovaSRadixSortDoubleALL::FrolovaSRadixSortDoubleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleALL::ValidationImpl() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int res = 0;
  if (rank == 0) {
    res = GetInput().empty() ? 0 : 1;
  }
  MPI_Bcast(&res, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return res == 1;
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
  MPI_Scatterv(rank == 0 ? GetInput().data() : nullptr, sendcounts.data(), displs.data(), MPI_DOUBLE, local_data.data(),
               sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  LocalRadixSort(local_data);

  std::vector<double> gathered;
  if (rank == 0) {
    gathered.resize(total_size);
  }

  MPI_Gatherv(local_data.data(), sendcounts[rank], MPI_DOUBLE, gathered.data(), sendcounts.data(), displs.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
