#include "lobanov_d_multi_matrix_crs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <thread>
#include <utility>
#include <vector>

#include "util/include/util.hpp"

namespace lobanov_d_multi_matrix_crs {

namespace {
constexpr double kEpsilon = 1e-15;
}  // namespace

LobanovMultyMatrixSTL::LobanovMultyMatrixSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LobanovMultyMatrixSTL::ValidationImpl() {
  const auto &input_data = GetInput();
  const auto &a = input_data.first;
  const auto &b = input_data.second;
  return a.column_count == b.row_count && a.row_count > 0 && b.column_count > 0;
}

bool LobanovMultyMatrixSTL::PreProcessingImpl() {
  return true;
}

void LobanovMultyMatrixSTL::SortIndices(std::vector<int> &vec) {
  std::ranges::sort(vec);  // Используйте ranges::sort
}

CompressedRowMatrix LobanovMultyMatrixSTL::MultiplyMatrices(const CompressedRowMatrix &a,
                                                            const CompressedRowMatrix &b) {
  CompressedRowMatrix result;
  result.row_count = a.row_count;
  result.column_count = b.column_count;
  result.row_pointer_data.assign(a.row_count + 1, 0);
  
  // Транспонирование b
  CompressedRowMatrix bt;
  bt.row_count = b.column_count;
  bt.column_count = b.row_count;
  bt.row_pointer_data.assign(bt.row_count + 1, 0);
  
  for (int col : b.column_index_data) {
    bt.row_pointer_data[col + 1]++;
  }
  for (int i = 0; i < bt.row_count; ++i) {
    bt.row_pointer_data[i + 1] += bt.row_pointer_data[i];
  }
  
  bt.value_data.resize(b.value_data.size());
  bt.column_index_data.resize(b.column_index_data.size());
  
  std::vector<int> current_pos = bt.row_pointer_data;
  for (int i = 0; i < b.row_count; ++i) {
    for (int j = b.row_pointer_data[i]; j < b.row_pointer_data[i + 1]; ++j) {
      int col = b.column_index_data[j];
      int dest = current_pos[col]++;
      bt.value_data[dest] = b.value_data[j];
      bt.column_index_data[dest] = i;
    }
  }
  
  // Временные хранилища для результата
  std::vector<std::vector<double>> temp_values(a.row_count);
  std::vector<std::vector<int>> temp_cols(a.row_count);
  
  int num_threads = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;
  
  auto worker = [&](int start_row, int end_row) {
    for (int i = start_row; i < end_row; ++i) {
      std::vector<double> vals;
      std::vector<int> cols;
      
      for (int j = 0; j < bt.row_count; ++j) {
        double sum = 0.0;
        int pa = a.row_pointer_data[i];
        int pb = bt.row_pointer_data[j];
        const int ea = a.row_pointer_data[i + 1];
        const int eb = bt.row_pointer_data[j + 1];
        
        while (pa < ea && pb < eb) {
          if (a.column_index_data[pa] == bt.column_index_data[pb]) {
            sum += a.value_data[pa] * bt.value_data[pb];
            ++pa;
            ++pb;
          } else if (a.column_index_data[pa] < bt.column_index_data[pb]) {
            ++pa;
          } else {
            ++pb;
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
  
  int rows_per_thread = a.row_count / num_threads;
  int remainder = a.row_count % num_threads;
  int current_start = 0;
  
  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    int current_end = current_start + rows_per_thread + (thread_idx < remainder ? 1 : 0);
    if (current_start < current_end) {
      threads.emplace_back(worker, current_start, current_end);
    }
    current_start = current_end;
  }
  
  for (auto &thread : threads) {
    thread.join();
  }
  
  // Сборка результата
  for (int i = 0; i < a.row_count; ++i) {
    result.row_pointer_data[i + 1] = result.row_pointer_data[i] + static_cast<int>(temp_values[i].size());
  }
  
  int total_nz = result.row_pointer_data[a.row_count];
  result.value_data.reserve(total_nz);
  result.column_index_data.reserve(total_nz);
  result.non_zero_count = total_nz;
  
  for (int i = 0; i < a.row_count; ++i) {
    result.value_data.insert(result.value_data.end(), temp_values[i].begin(), temp_values[i].end());
    result.column_index_data.insert(result.column_index_data.end(), temp_cols[i].begin(), temp_cols[i].end());
  }
  
  return result;
}

bool LobanovMultyMatrixSTL::RunImpl() {
  const auto &input_data = GetInput();
  const auto &a = input_data.first;
  const auto &b = input_data.second;
  
  CompressedRowMatrix c = MultiplyMatrices(a, b);
  GetOutput() = c;
  
  return true;
}

bool LobanovMultyMatrixSTL::PostProcessingImpl() {
  return true;
}

}  // namespace lobanov_d_multi_matrix_crs