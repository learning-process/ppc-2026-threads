#include "liulin_y_complex_ccs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <map>
#include <thread>
#include <utility>
#include <vector>

#include "liulin_y_complex_ccs/common/include/common.hpp"

namespace liulin_y_complex_ccs {

namespace {

constexpr double kEpsilon = 1e-10;

using LocalDict = std::map<std::pair<int, int>, std::complex<double>>;

bool IsValidCCS(const CCSMatrix &mat) {
  if (mat.count_rows <= 0 || mat.count_cols <= 0) {
    return false;
  }
  if (mat.col_index.size() != static_cast<size_t>(mat.count_cols) + 1) {
    return false;
  }
  if (mat.col_index[0] != 0) {
    return false;
  }
  if (mat.values.size() != mat.row_index.size()) {
    return false;
  }
  if (std::cmp_not_equal(mat.col_index.back(), mat.values.size())) {
    return false;
  }
  return true;
}

void TransposeCCS(const CCSMatrix &mat, CCSMatrix &mat_t) {
  mat_t.count_rows = mat.count_cols;
  mat_t.count_cols = mat.count_rows;
  const size_t nnz = mat.values.size();
  mat_t.values.resize(nnz);
  mat_t.row_index.resize(nnz);
  mat_t.col_index.assign(static_cast<size_t>(mat_t.count_cols) + 1, 0);

  for (size_t i = 0; i < nnz; ++i) {
    mat_t.col_index[static_cast<size_t>(mat.row_index[i]) + 1]++;
  }
  for (int i = 0; i < mat_t.count_cols; ++i) {
    mat_t.col_index[static_cast<size_t>(i) + 1] += mat_t.col_index[static_cast<size_t>(i)];
  }

  std::vector<int> current_pos(mat_t.col_index.begin(), mat_t.col_index.end());
  for (int j = 0; j < mat.count_cols; ++j) {
    for (int k = mat.col_index[static_cast<size_t>(j)]; k < mat.col_index[static_cast<size_t>(j) + 1]; ++k) {
      const int row = mat.row_index[static_cast<size_t>(k)];
      const int dest = current_pos[static_cast<size_t>(row)]++;
      mat_t.row_index[static_cast<size_t>(dest)] = j;
      mat_t.values[static_cast<size_t>(dest)] = mat.values[static_cast<size_t>(k)];
    }
  }
}

std::vector<int> BuildColumnIndices(const CCSMatrix &mat_a) {
  std::vector<int> a_cols(mat_a.values.size());
  for (int j = 0; j < mat_a.count_cols; ++j) {
    for (int k = mat_a.col_index[static_cast<size_t>(j)]; k < mat_a.col_index[static_cast<size_t>(j) + 1]; ++k) {
      a_cols[static_cast<size_t>(k)] = j;
    }
  }
  return a_cols;
}

void ProcessChunk(int start_idx, int end_idx, const CCSMatrix &mat_a, const std::vector<int> &a_cols,
                  const CCSMatrix &mat_b_t, LocalDict &local_map) {
  for (int k = start_idx; k < end_idx; ++k) {
    const int row_left = mat_a.row_index[static_cast<size_t>(k)];
    const int col_left = a_cols[static_cast<size_t>(k)];
    const std::complex<double> val_left = mat_a.values[static_cast<size_t>(k)];

    const int b_start = mat_b_t.col_index[static_cast<size_t>(col_left)];
    const int b_end = mat_b_t.col_index[static_cast<size_t>(col_left) + 1];

    for (int idx_p = b_start; idx_p < b_end; ++idx_p) {
      const int col_right = mat_b_t.row_index[static_cast<size_t>(idx_p)];
      const std::complex<double> val_right = mat_b_t.values[static_cast<size_t>(idx_p)];

      local_map[{col_right, row_left}] += val_left * val_right;
    }
  }
}

void AggregateLocalMaps(const std::vector<LocalDict> &local_maps, LocalDict &global_map) {
  for (const auto &local_map : local_maps) {
    for (const auto &[key, value] : local_map) {
      global_map[key] += value;
    }
  }
}

void ConstructResultMatrix(const LocalDict &global_map, int rows, int cols, CCSMatrix &mat_res) {
  mat_res.count_rows = rows;
  mat_res.count_cols = cols;
  mat_res.values.clear();
  mat_res.row_index.clear();
  mat_res.col_index.assign(static_cast<size_t>(cols) + 1, 0);

  for (const auto &[key, value] : global_map) {
    if (std::abs(value.real()) > kEpsilon || std::abs(value.imag()) > kEpsilon) {
      const int col = key.first;
      const int row = key.second;

      mat_res.values.push_back(value);
      mat_res.row_index.push_back(row);
      mat_res.col_index[static_cast<size_t>(col) + 1]++;
    }
  }

  for (int i = 0; i < cols; ++i) {
    mat_res.col_index[static_cast<size_t>(i) + 1] += mat_res.col_index[static_cast<size_t>(i)];
  }
}

}  // namespace

LiulinYComplexCcsStl::LiulinYComplexCcsStl(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LiulinYComplexCcsStl::ValidationImpl() {
  const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;
  return IsValidCCS(mat_a) && IsValidCCS(mat_b) && (mat_a.count_cols == mat_b.count_rows);
}

bool LiulinYComplexCcsStl::PreProcessingImpl() {
  return true;
}

bool LiulinYComplexCcsStl::RunImpl() {
  const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;
  auto &mat_res = GetOutput();

  const int nnz_a = static_cast<int>(mat_a.values.size());
  if (nnz_a == 0) {
    mat_res.count_rows = mat_a.count_rows;
    mat_res.count_cols = mat_b.count_cols;
    mat_res.col_index.assign(static_cast<size_t>(mat_res.count_cols) + 1, 0);
    return true;
  }

  CCSMatrix mat_b_t;
  TransposeCCS(mat_b, mat_b_t);
  std::vector<int> a_cols = BuildColumnIndices(mat_a);

  unsigned int hw = std::thread::hardware_concurrency();
  int num_threads = std::min(nnz_a, (hw == 0) ? 1 : static_cast<int>(hw));

  std::vector<std::thread> threads;
  std::vector<LocalDict> local_maps(static_cast<size_t>(num_threads));
  int chunk = nnz_a / num_threads;
  int rem = nnz_a % num_threads;

  int start = 0;
  for (int i = 0; i < num_threads; ++i) {
    int end = start + chunk + (i < rem ? 1 : 0);
    threads.emplace_back(ProcessChunk, start, end, std::ref(mat_a), std::ref(a_cols), std::ref(mat_b_t),
                         std::ref(local_maps[static_cast<size_t>(i)]));
    start = end;
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  LocalDict global_map;
  AggregateLocalMaps(local_maps, global_map);
  ConstructResultMatrix(global_map, mat_a.count_rows, mat_b.count_cols, mat_res);

  return true;
}

bool LiulinYComplexCcsStl::PostProcessingImpl() {
  return true;
}

}  // namespace liulin_y_complex_ccs
