#include "popova_e_radix_sort_for_double_with_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <random>
#include <utility>
#include <vector>

#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

namespace {

uint64_t DoubleToSortable(double value) {
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(double));
  bool is_negative = (bits >> 63) == 1;
  if (is_negative) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);
  }
  return bits;
}

double SortableToDouble(uint64_t bits) {
  bool is_negative = (bits >> 63) == 1;
  if (is_negative) {
    bits ^= (1ULL << 63);
  } else {
    bits = ~bits;
  }
  double value = 0;
  std::memcpy(&value, &bits, sizeof(double));
  return value;
}

void RadixSortUInt(std::vector<uint64_t> &arr) {
  if (arr.empty()) {
    return;
  }

  const int bytes_count = 8;
  const int base = 256;
  std::vector<uint64_t> buffer(arr.size());

  for (int byte_index = 0; byte_index < bytes_count; byte_index++) {
    int sdvig = byte_index * 8;
    std::array<size_t, base> count = {0};

    for (const auto &val : arr) {
      count.at((val >> sdvig) & 0xFF)++;
    }

    size_t offset = 0;
    for (auto &c : count) {
      size_t tmp = c;
      c = offset;
      offset += tmp;
    }

    for (const auto &val : arr) {
      size_t pos = (val >> sdvig) & 0xFF;
      buffer.at(count.at(pos)) = val;
      count.at(pos)++;
    }
    arr = buffer;
  }
}

std::vector<double> MergeSorted(const std::vector<double> &left, const std::vector<double> &right) {
  std::vector<double> res;
  res.reserve(left.size() + right.size());
  size_t i = 0;
  size_t j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      res.push_back(left[i++]);
    } else {
      res.push_back(right[j++]);
    }
  }
  while (i < left.size()) {
    res.push_back(left[i++]);
  }
  while (j < right.size()) {
    res.push_back(right[j++]);
  }
  return res;
}

double RandomDouble(double min_val, double max_val) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(min_val, max_val);
  return dis(gen);
}

bool IsSorted(const std::vector<double> &arr) {
  for (size_t i = 1; i < arr.size(); i++) {
    if (arr[i - 1] > arr[i]) {
      return false;
    }
  }
  return true;
}

bool SameData(const std::vector<double> &original, const std::vector<double> &result) {
  uint64_t hash_original = 0;
  uint64_t hash_result = 0;

  for (const double &value : original) {
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(double));
    hash_original ^= bits;
  }

  for (const double &value : result) {
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(double));
    hash_result ^= bits;
  }

  return hash_original == hash_result;
}

void SplitMpiData(int total_elements, int size, std::vector<int> &elem_count, std::vector<int> &start_index) {
  int elems_per_proc = total_elements / size;
  int ostat = total_elements % size;
  int curr_index = 0;

  for (int i = 0; i < size; i++) {
    if (i < ostat) {
      elem_count[i] = elems_per_proc + 1;
    } else {
      elem_count[i] = elems_per_proc;
    }
    start_index[i] = curr_index;
    curr_index += elem_count[i];
  }
}

void MergeChunks(const std::vector<double> &gathered, const std::vector<int> &counts, const std::vector<int> &starts,
                 int size, std::vector<double> &output) {
  output.clear();
  for (int i = 0; i < size; i++) {
    std::vector<double> piece(gathered.begin() + starts[i], gathered.begin() + starts[i] + counts[i]);
    if (!piece.empty()) {
      if (output.empty()) {
        output = std::move(piece);
      } else {
        output = MergeSorted(output, piece);
      }
    }
  }
}

}  // namespace

PopovaERadixSorForDoubleWithSimpleMergeALL::PopovaERadixSorForDoubleWithSimpleMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool PopovaERadixSorForDoubleWithSimpleMergeALL::ValidationImpl() {
  return GetInput() > 0;
}

bool PopovaERadixSorForDoubleWithSimpleMergeALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    int size = GetInput();
    array_.resize(size);
    for (int i = 0; i < size; i++) {
      array_[i] = RandomDouble(-100.0, 100.0);
    }
  }
  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int total_elements = 0;
  if (rank == 0) {
    total_elements = static_cast<int>(array_.size());
  }
  MPI_Bcast(&total_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> elem_count(size);
  std::vector<int> start_index(size);
  SplitMpiData(total_elements, size, elem_count, start_index);

  int my_elem_count = elem_count[rank];
  std::vector<double> local_array(my_elem_count);

  MPI_Scatterv(array_.data(), elem_count.data(), start_index.data(), MPI_DOUBLE, local_array.data(), my_elem_count,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::vector<uint64_t> local_bits(my_elem_count);

#pragma omp parallel default(none) shared(my_elem_count, local_bits, local_array)
  {
#pragma omp for
    for (int i = 0; i < my_elem_count; i++) {
      local_bits[i] = DoubleToSortable(local_array[i]);
    }
  }

  RadixSortUInt(local_bits);

#pragma omp parallel default(none) shared(my_elem_count, local_array, local_bits)
  {
#pragma omp for
    for (int i = 0; i < my_elem_count; i++) {
      local_array[i] = SortableToDouble(local_bits[i]);
    }
  }

  std::vector<double> gathered_array;
  if (rank == 0) {
    gathered_array.resize(total_elements);
  }

  MPI_Gatherv(local_array.data(), my_elem_count, MPI_DOUBLE, gathered_array.data(), elem_count.data(),
              start_index.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    MergeChunks(gathered_array, elem_count, start_index, size, result_);
  }

  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int is_ok = 0;
  if (rank == 0) {
    bool sorted = IsSorted(result_);
    bool same = SameData(array_, result_);
    if (sorted && same) {
      is_ok = 1;
    } else {
      is_ok = 0;
    }
    GetOutput() = is_ok;
  }

  MPI_Bcast(&is_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    GetOutput() = is_ok;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads
