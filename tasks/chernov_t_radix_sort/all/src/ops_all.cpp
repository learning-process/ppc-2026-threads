#include "chernov_t_radix_sort/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"

namespace chernov_t_radix_sort {

// ============================================================================
// Конструктор и базовые методы
// ============================================================================

ChernovTRadixSortALL::ChernovTRadixSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ChernovTRadixSortALL::ValidationImpl() {
  return true;
}

bool ChernovTRadixSortALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

// ============================================================================
// Константы
// ============================================================================

constexpr int kBitsPerDigit = 8;
constexpr int kRadix = 1 << kBitsPerDigit;  // 256
constexpr uint32_t kSignMask = 0x80000000U;

// ============================================================================
// Локальная сортировка с OpenMP (без гонок данных)
// ============================================================================

void ChernovTRadixSortALL::RadixSortLSDParallelOMP(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }
  const size_t n = data.size();

  // 1. Преобразование signed -> unsigned с инверсией знака
  std::vector<uint32_t> temp(n);
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    temp[i] = static_cast<uint32_t>(data[i]) ^ kSignMask;
  }

  std::vector<uint32_t> buffer(n);
  int num_threads = omp_get_max_threads();

  // 2. 4 прохода по байтам
  for (int byte = 0; byte < 4; ++byte) {
    const int shift = byte * kBitsPerDigit;

    // Локальные гистограммы для каждого потока
    std::vector<std::vector<int>> local_counts(static_cast<size_t>(num_threads), std::vector<int>(kRadix, 0));

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; ++i) {
      int tid = omp_get_thread_num();
      int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
      local_counts[static_cast<size_t>(tid)][static_cast<size_t>(digit)]++;
    }

    // Суммируем локальные гистограммы в глобальную
    std::vector<int> global_count(kRadix, 0);
    for (int t = 0; t < num_threads; ++t) {
      for (int d = 0; d < kRadix; ++d) {
        global_count[d] += local_counts[static_cast<size_t>(t)][static_cast<size_t>(d)];
      }
    }

    // Префиксные суммы (определяем начальные позиции)
    for (int i = 1; i < kRadix; ++i) {
      global_count[i] += global_count[i - 1];
    }

    // ================================================================
    // РАССЕИВАНИЕ — ПОСЛЕДОВАТЕЛЬНОЕ (без гонок данных)
    // ================================================================
    std::vector<int> count = global_count;  // копия для модификации
    for (size_t i = n; i-- > 0;) {
      int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
      buffer[static_cast<size_t>(--count[static_cast<size_t>(digit)])] = temp[i];
    }

    temp.swap(buffer);
  }

  // 3. Обратное преобразование
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    data[i] = static_cast<int>(temp[i] ^ kSignMask);
  }
}

// ============================================================================
// Простое слияние двух отсортированных массивов
// ============================================================================

void ChernovTRadixSortALL::SimpleMerge(const std::vector<int> &left, const std::vector<int> &right,
                                       std::vector<int> &result) {
  result.resize(left.size() + right.size());
  std::merge(left.begin(), left.end(), right.begin(), right.end(), result.begin());
}

// ============================================================================
// RunImpl: MPI + OMP (гибрид)
// ============================================================================

bool ChernovTRadixSortALL::RunImpl() {
  auto &input_data = GetInput();
  int rank = -1, size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const size_t total_elements = input_data.size();

  // Пустой массив — сразу выходим
  if (total_elements == 0) {
    GetOutput().clear();
    return true;
  }

  // 1. Расчёт размеров чанков (только на rank 0)
  std::vector<int> recv_counts(size, 0);
  std::vector<int> displs(size, 0);

  if (rank == 0) {
    int base = static_cast<int>(total_elements / static_cast<size_t>(size));
    int remainder = static_cast<int>(total_elements % static_cast<size_t>(size));
    int current_disp = 0;
    for (int i = 0; i < size; ++i) {
      recv_counts[i] = base + (i < remainder ? 1 : 0);
      displs[i] = current_disp;
      current_disp += recv_counts[i];
    }
  }

  // 2. Рассылаем размеры чанков всем процессам
  int local_n = 0;
  MPI_Scatter(recv_counts.data(), 1, MPI_INT, &local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // 3. Рассылаем данные
  std::vector<int> local_data;
  if (local_n > 0) {
    local_data.resize(static_cast<size_t>(local_n));
    MPI_Scatterv(input_data.data(), recv_counts.data(), displs.data(), MPI_INT, local_data.data(), local_n, MPI_INT, 0,
                 MPI_COMM_WORLD);

    // 4. Локальная сортировка с OpenMP
    RadixSortLSDParallelOMP(local_data);
  }

  // 5. Собираем результаты на rank 0
  std::vector<int> global_result;
  if (rank == 0) {
    global_result.resize(total_elements);
  }

  MPI_Gatherv(local_data.empty() ? nullptr : local_data.data(), local_n, MPI_INT, global_result.data(),
              recv_counts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  // 6. Слияние всех чанков на rank 0
  if (rank == 0) {
    std::vector<int> merged;
    int offset = 0;

    // Находим первый непустой чанк
    int start_idx = -1;
    for (int i = 0; i < size; ++i) {
      if (recv_counts[i] > 0) {
        start_idx = i;
        offset = displs[i];
        break;
      }
    }

    if (start_idx != -1) {
      // Инициализируем merged первым непустым чанком
      merged.assign(global_result.begin() + offset, global_result.begin() + offset + recv_counts[start_idx]);
      offset += recv_counts[start_idx];

      // Последовательно сливаем остальные чанки
      for (int i = start_idx + 1; i < size; ++i) {
        if (recv_counts[i] > 0) {
          std::vector<int> next_part(global_result.begin() + offset, global_result.begin() + offset + recv_counts[i]);
          std::vector<int> new_merged;
          SimpleMerge(merged, next_part, new_merged);
          merged = std::move(new_merged);
          offset += recv_counts[i];
        }
      }
    }
    GetOutput() = std::move(merged);
  } else {
    GetOutput().clear();
  }

  // 7. Рассылаем итоговый результат всем процессам
  //    (чтобы PostProcessingImpl работал на всех рангах)
  int out_size = static_cast<int>(GetOutput().size());
  MPI_Bcast(&out_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    GetOutput().resize(static_cast<size_t>(out_size));
  }

  if (out_size > 0) {
    MPI_Bcast(GetOutput().data(), out_size, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return true;
}

// ============================================================================
// PostProcessing: проверка, что результат отсортирован
// ============================================================================

bool ChernovTRadixSortALL::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace chernov_t_radix_sort
