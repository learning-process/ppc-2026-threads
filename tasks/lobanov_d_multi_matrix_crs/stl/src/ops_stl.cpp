#include "lobanov_d_multi_matrix_crs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

#include "util/include/util.hpp"

namespace lobanov_d_multi_matrix_crs {

namespace {
constexpr double kEpsilon = 1e-15;
}

LobanovMultyMatrixSTL::LobanovMultyMatrixSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LobanovMultyMatrixSTL::ValidationImpl() {
  const auto &input_data = GetInput();
  const auto &A = input_data.first;
  const auto &B = input_data.second;
  return A.column_count == B.row_count && A.row_count > 0 && B.column_count > 0;
}

bool LobanovMultyMatrixSTL::PreProcessingImpl() {
  return true;
}

void LobanovMultyMatrixSTL::SortIndices(std::vector<int> &vec) {
  std::sort(vec.begin(), vec.end());
}

CompressedRowMatrix LobanovMultyMatrixSTL::MultiplyMatrices(const CompressedRowMatrix &A,
                                                            const CompressedRowMatrix &B) {
  CompressedRowMatrix result;
  result.row_count = A.row_count;
  result.column_count = B.column_count;
  result.row_pointer_data.assign(A.row_count + 1, 0);

  // Транспонирование B
  CompressedRowMatrix Bt;
  Bt.row_count = B.column_count;
  Bt.column_count = B.row_count;
  Bt.row_pointer_data.assign(Bt.row_count + 1, 0);

  for (int col : B.column_index_data) {
    Bt.row_pointer_data[col + 1]++;
  }
  for (int i = 0; i < Bt.row_count; ++i) {
    Bt.row_pointer_data[i + 1] += Bt.row_pointer_data[i];
  }

  Bt.value_data.resize(B.value_data.size());
  Bt.column_index_data.resize(B.column_index_data.size());

  std::vector<int> current_pos = Bt.row_pointer_data;
  for (int i = 0; i < B.row_count; ++i) {
    for (int j = B.row_pointer_data[i]; j < B.row_pointer_data[i + 1]; ++j) {
      int col = B.column_index_data[j];
      int dest = current_pos[col]++;
      Bt.value_data[dest] = B.value_data[j];
      Bt.column_index_data[dest] = i;
    }
  }

  // Временные хранилища для результата
  std::vector<std::vector<double>> temp_values(A.row_count);
  std::vector<std::vector<int>> temp_cols(A.row_count);

  int num_threads = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;

  auto worker = [&](int start_row, int end_row) {
    for (int i = start_row; i < end_row; ++i) {
      std::vector<double> vals;
      std::vector<int> cols;

      for (int j = 0; j < Bt.row_count; ++j) {
        double sum = 0.0;
        int pa = A.row_pointer_data[i];
        int pb = Bt.row_pointer_data[j];
        const int ea = A.row_pointer_data[i + 1];
        const int eb = Bt.row_pointer_data[j + 1];

        while (pa < ea && pb < eb) {
          if (A.column_index_data[pa] == Bt.column_index_data[pb]) {
            sum += A.value_data[pa] * Bt.value_data[pb];
            pa++;
            pb++;
          } else if (A.column_index_data[pa] < Bt.column_index_data[pb]) {
            pa++;
          } else {
            pb++;
          }
        }

        if (std::abs(sum) > kEpsilon) {
          vals.push_back(sum);
          cols.push_back(j);
        }
      }
      temp_values[i] = std::move(vals);
      temp_cols[i] = std::move(cols);
    }
  };

  int rows_per_thread = A.row_count / num_threads;
  int remainder = A.row_count % num_threads;
  int current_start = 0;

  for (int t = 0; t < num_threads; ++t) {
    int current_end = current_start + rows_per_thread + (t < remainder ? 1 : 0);
    if (current_start < current_end) {
      threads.emplace_back(worker, current_start, current_end);
    }
    current_start = current_end;
  }

  for (auto &thread : threads) {
    thread.join();
  }

  // Сборка результата
  for (int i = 0; i < A.row_count; ++i) {
    result.row_pointer_data[i + 1] = result.row_pointer_data[i] + static_cast<int>(temp_values[i].size());
  }

  int total_nz = result.row_pointer_data[A.row_count];
  result.value_data.reserve(total_nz);
  result.column_index_data.reserve(total_nz);
  result.non_zero_count = total_nz;

  for (int i = 0; i < A.row_count; ++i) {
    result.value_data.insert(result.value_data.end(), temp_values[i].begin(), temp_values[i].end());
    result.column_index_data.insert(result.column_index_data.end(), temp_cols[i].begin(), temp_cols[i].end());
  }

  return result;
}

bool LobanovMultyMatrixSTL::RunImpl() {
  const auto &input_data = GetInput();
  const auto &A = input_data.first;
  const auto &B = input_data.second;

  CompressedRowMatrix C = MultiplyMatrices(A, B);
  GetOutput() = C;

  return true;
}

bool LobanovMultyMatrixSTL::PostProcessingImpl() {
  return true;
}

}  // namespace lobanov_d_multi_matrix_crs
