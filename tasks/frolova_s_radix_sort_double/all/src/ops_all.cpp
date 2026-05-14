#include "frolova_s_radix_sort_double/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"

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

  // Основной цикл сортировки по байтам
  for (int pass = 0; pass < num_passes; pass++) {
    std::vector<int> count(radix, 0);
    for (double value : chunk) {
      auto bits = std::bit_cast<uint64_t>(value);
      int byte = static_cast<int>((bits >> (pass * num_bits)) & 0xFF);
      count[byte]++;
    }

    int total = 0;
    for (int i = 0; i < radix; i++) {
      int old = count[i];
      count[i] = total;
      total += old;
    }

    for (double value : chunk) {
      auto bits = std::bit_cast<uint64_t>(value);
      int byte = static_cast<int>((bits >> (pass * num_bits)) & 0xFF);
      temp[count[byte]++] = value;
    }
    chunk.swap(temp);
  }

  // Разделение на отрицательные и положительные числа для исправления порядка
  std::vector<double> negative;
  std::vector<double> positive;
  negative.reserve(chunk.size());
  positive.reserve(chunk.size());

  for (double val : chunk) {
    if ((std::bit_cast<uint64_t>(val) >> 63) != 0U) {
      negative.push_back(val);
    } else {
      positive.push_back(val);
    }
  }

  // Отрицательные числа хранятся в обратном порядке при поразрядной сортировке битов IEEE 754
  std::ranges::reverse(negative);

  size_t pos = 0;
  for (double val : negative) {
    chunk[pos++] = val;
  }
  for (double val : positive) {
    chunk[pos++] = val;
  }
}

}  // namespace

FrolovaSRadixSortDoubleALL::FrolovaSRadixSortDoubleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return !GetInput().empty();
  }
  return true;
}

bool FrolovaSRadixSortDoubleALL::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int total_size = 0;
  if (rank == 0) {
    total_size = static_cast<int>(GetInput().size());
  }

  // Рассылаем общий размер массива всем процессам
  MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total_size == 0) {
    return false;
  }

  // Вычисляем размеры блоков для каждого процесса
  std::vector<int> sendcounts(size);
  std::vector<int> displs(size);
  int remainder = total_size % size;
  int offset = 0;

  for (int i = 0; i < size; ++i) {
    sendcounts[i] = total_size / size + (i < remainder ? 1 : 0);
    displs[i] = offset;
    offset += sendcounts[i];
  }

  // Выделяем память под локальный блок
  std::vector<double> local_data(sendcounts[rank]);

  // Распределяем данные
  if (rank == 0) {
    MPI_Scatterv(GetInput().data(), sendcounts.data(), displs.data(), MPI_DOUBLE, local_data.data(), sendcounts[rank],
                 MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else {
    MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE, local_data.data(), sendcounts[rank], MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
  }

  // Локальная поразрядная сортировка каждого блока
  LocalRadixSort(local_data);

  // Сбор данных обратно на 0-й процесс
  std::vector<double> gathered_data;
  if (rank == 0) {
    gathered_data.resize(total_size);
    MPI_Gatherv(local_data.data(), sendcounts[rank], MPI_DOUBLE, gathered_data.data(), sendcounts.data(), displs.data(),
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else {
    MPI_Gatherv(local_data.data(), sendcounts[rank], MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
  }

  // Простое слияние на главном процессе (rank 0)
  if (rank == 0) {
    std::vector<double> merged_result;
    merged_result.assign(gathered_data.begin(), gathered_data.begin() + sendcounts[0]);

    for (int i = 1; i < size; ++i) {
      if (sendcounts[i] == 0) {
        continue;
      }
      std::vector<double> merged(merged_result.size() + sendcounts[i]);
      auto next_chunk_begin = gathered_data.begin() + displs[i];
      auto next_chunk_end = next_chunk_begin + sendcounts[i];

      std::merge(merged_result.begin(), merged_result.end(), next_chunk_begin, next_chunk_end, merged.begin());

      merged_result = std::move(merged);
    }

    GetOutput() = std::move(merged_result);
  }

  return true;
}

bool FrolovaSRadixSortDoubleALL::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
