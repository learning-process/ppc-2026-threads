#include "romanov_m_matrix_ccs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "romanov_m_matrix_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace romanov_m_matrix_ccs {

RomanovMMatrixCCSSTL::RomanovMMatrixCCSSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RomanovMMatrixCCSSTL::ValidationImpl() {
  return GetInput().first.cols_num == GetInput().second.rows_num;
}

bool RomanovMMatrixCCSSTL::PreProcessingImpl() {
  return true;
}

void RomanovMMatrixCCSSTL::MultiplyColumn(size_t col_index, const MatrixCCS &a, const MatrixCCS &b,
                                          std::vector<double> &temp_v, std::vector<size_t> &temp_r) {
  std::vector<double> accumulator(a.rows_num, 0.0);
  std::vector<bool> row_mask(a.rows_num, false);
  std::vector<size_t> active_rows;

  for (size_t kb = b.col_ptrs[col_index]; kb < b.col_ptrs[col_index + 1]; ++kb) {
    size_t k = b.row_inds[kb];
    double v_b = b.vals[kb];

    for (size_t ka = a.col_ptrs[k]; ka < a.col_ptrs[k + 1]; ++ka) {
      size_t i = a.row_inds[ka];
      if (!row_mask[i]) {
        row_mask[i] = true;
        active_rows.push_back(i);
      }
      accumulator[i] += a.vals[ka] * v_b;
    }
  }

  std::ranges::sort(active_rows);
  for (size_t row_idx : active_rows) {
    if (std::abs(accumulator[row_idx]) > 1e-12) {
      temp_v.push_back(accumulator[row_idx]);
      temp_r.push_back(row_idx);
    }
  }
}

bool RomanovMMatrixCCSSTL::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput();

  c.rows_num = a.rows_num;
  c.cols_num = b.cols_num;
  c.col_ptrs.assign(c.cols_num + 1, 0);

  std::vector<std::vector<double>> temp_vals(c.cols_num);
  std::vector<std::vector<size_t>> temp_rows(c.cols_num);

  unsigned int num_threads = ppc::util::GetNumThreads();
  if (num_threads == 0) {
    num_threads = 1;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  auto worker = [&](size_t start, size_t end) {
    for (size_t j = start; j < end; ++j) {
      MultiplyColumn(j, a, b, temp_vals[j], temp_rows[j]);
    }
  };

  size_t chunk_size = b.cols_num / num_threads;
  size_t remainder = b.cols_num % num_threads;
  size_t current_start = 0;

  for (unsigned int i = 0; i < num_threads; ++i) {
    size_t current_end = current_start + chunk_size + (i < remainder ? 1 : 0);
    if (current_start < current_end) {
      threads.emplace_back(worker, current_start, current_end);
    }
    current_start = current_end;
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  size_t total_nnz = 0;
  for (size_t j = 0; j < b.cols_num; ++j) {
    c.col_ptrs[j] = total_nnz;
    total_nnz += temp_vals[j].size();
  }
  c.col_ptrs[b.cols_num] = total_nnz;
  c.nnz = total_nnz;

  c.vals.reserve(total_nnz);
  c.row_inds.reserve(total_nnz);
  for (size_t j = 0; j < b.cols_num; ++j) {
    c.vals.insert(c.vals.end(), temp_vals[j].begin(), temp_vals[j].end());
    c.row_inds.insert(c.row_inds.end(), temp_rows[j].begin(), temp_rows[j].end());
  }

  return true;
}

bool RomanovMMatrixCCSSTL::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_m_matrix_ccs
