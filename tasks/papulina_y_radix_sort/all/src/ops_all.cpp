#include "papulina_y_radix_sort/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"

namespace papulina_y_radix_sort {

PapulinaYRadixSortALL::PapulinaYRadixSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool PapulinaYRadixSortALL::ValidationImpl() {
  return true;
}

bool PapulinaYRadixSortALL::PreProcessingImpl() {
  return true;
}
std::vector<double> PapulinaYRadixSortALL::SimpleMerge(const std::vector<double> &a, const std::vector<double> &b) {
  std::vector<double> res;
  res.reserve(a.size() + b.size());
  size_t i = 0;
  size_t j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      res.push_back(a[i++]);
    } else {
      res.push_back(b[j++]);
    }
  }
  while (i < a.size()) {
    res.push_back(a[i++]);
  }
  while (j < b.size()) {
    res.push_back(b[j++]);
  }
  return res;
}

bool PapulinaYRadixSortALL::RunImpl() {
  int rank = 0;
  int size_comm = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size_comm);

  int total_size = 0;
  if (rank == 0) {
    total_size = static_cast<int>(GetInput().size());
  }
  MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total_size == 0) {
    return true;
  }

  std::vector<int> send_counts(size_comm);
  std::vector<int> offsets(size_comm);
  int chunk = total_size / size_comm;
  int remainder = total_size % size_comm;

  for (int i = 0; i < size_comm; ++i) {
    send_counts[i] = chunk + (i < remainder ? 1 : 0);
    offsets[i] = (i == 0) ? 0 : offsets[i - 1] + send_counts[i - 1];
  }

  std::vector<double> local_data(send_counts[rank]);
  double *in_ptr = (rank == 0) ? GetInput().data() : nullptr;
  MPI_Scatterv(in_ptr, send_counts.data(), offsets.data(), MPI_DOUBLE, local_data.data(), send_counts[rank], MPI_DOUBLE,
               0, MPI_COMM_WORLD);

  RadixSort(local_data.data(), send_counts[rank]);

  if (rank == 0) {
    std::vector<double> final_res = local_data;
    for (int i = 1; i < size_comm; ++i) {
      if (send_counts[i] > 0) {
        std::vector<double> recv_buf(send_counts[i]);
        MPI_Recv(recv_buf.data(), send_counts[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        final_res = SimpleMerge(final_res, recv_buf);
      }
    }
    GetOutput() = final_res;
  } else if (send_counts[rank] > 0) {
    MPI_Send(local_data.data(), send_counts[rank], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  GetOutput().resize(total_size);
  MPI_Bcast(GetOutput().data(), total_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return true;
}

bool PapulinaYRadixSortALL::PostProcessingImpl() {
  return true;
}
void PapulinaYRadixSortALL::SortByByte(uint64_t *bytes, uint64_t *out, int byte, int size) {
  auto *byte_view = reinterpret_cast<unsigned char *>(bytes);
  std::array<int, 256> counter = {0};

  for (int i = 0; i < size; i++) {
    int index = byte_view[(8 * i) + byte];
    *(counter.data() + index) += 1;
  }

  int tmp = 0;
  for (int j = 0; j < 256; j++) {
    int a = *(counter.data() + j);
    *(counter.data() + j) = tmp;
    tmp += a;
  }

  for (int i = 0; i < size; i++) {
    int index = byte_view[(8 * i) + byte];
    out[*(counter.data() + index)] = bytes[i];
    *(counter.data() + index) += 1;
  }
}

uint64_t PapulinaYRadixSortALL::InBytes(double d) {
  uint64_t bits = 0;
  memcpy(&bits, &d, sizeof(double));
  if ((bits & kMask) != 0) {
    bits = ~bits;
  } else {
    bits = bits ^ kMask;
  }
  return bits;
}

double PapulinaYRadixSortALL::FromBytes(uint64_t bits) {
  double d = NAN;
  if ((bits & kMask) != 0) {
    bits = bits ^ kMask;
  } else {
    bits = ~bits;
  }
  memcpy(&d, &bits, sizeof(double));
  return d;
}

void PapulinaYRadixSortALL::RadixSort(double *arr, int size) {
  if (size <= 1) {
    return;
  }
  std::vector<uint64_t> bytes(size);
  std::vector<uint64_t> out(size);

#pragma omp parallel for default(none) shared(bytes, arr, size)
  for (int i = 0; i < size; i++) {
    bytes[i] = InBytes(arr[i]);
  }

  uint64_t *src_ptr = bytes.data();
  uint64_t *dst_ptr = out.data();

  for (int byte = 0; byte < 8; byte++) {
    SortByByte(src_ptr, dst_ptr, byte, size);
    std::swap(src_ptr, dst_ptr);
  }

#pragma omp parallel for default(none) shared(arr, src_ptr, size)
  for (int i = 0; i < size; i++) {
    arr[i] = FromBytes(src_ptr[i]);
  }
}
}  // namespace papulina_y_radix_sort
