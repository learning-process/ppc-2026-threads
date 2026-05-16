#include "potashnik_m_matrix_mult_complex/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <thread>
#include <utility>
#include <vector>

#include "potashnik_m_matrix_mult_complex/common/include/common.hpp"

namespace potashnik_m_matrix_mult_complex {

PotashnikMMatrixMultComplexSTL::PotashnikMMatrixMultComplexSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool PotashnikMMatrixMultComplexSTL::ValidationImpl() {
  const auto &matrix_left = std::get<0>(GetInput());
  const auto &matrix_right = std::get<1>(GetInput());
  return matrix_left.width == matrix_right.height;
}

bool PotashnikMMatrixMultComplexSTL::PreProcessingImpl() {
  return true;
}

namespace {

using Key = std::pair<size_t, size_t>;
using LocalMap = std::map<Key, Complex>;

void ProcessChunk(size_t begin, size_t end, const CCSMatrix &matrix_right, const std::vector<Complex> &val_left,
                  const std::vector<size_t> &row_ind_left, const std::vector<size_t> &col_ptr_left,
                  LocalMap &local_buffer) {
  for (size_t i = begin; i < end; ++i) {
    size_t row_left = row_ind_left[i];
    size_t col_left = col_ptr_left[i];
    Complex left_val = val_left[i];

    for (size_t j = 0; j < matrix_right.Count(); ++j) {
      if (col_left == matrix_right.row_ind[j]) {
        local_buffer[{row_left, matrix_right.col_ptr[j]}] += left_val * matrix_right.val[j];
      }
    }
  }
}

}  // namespace

bool PotashnikMMatrixMultComplexSTL::RunImpl() {
  const auto &matrix_left = std::get<0>(GetInput());
  const auto &matrix_right = std::get<1>(GetInput());

  const auto &val_left = matrix_left.val;
  const auto &row_ind_left = matrix_left.row_ind;
  const auto &col_ptr_left = matrix_left.col_ptr;
  size_t height_left = matrix_left.height;
  size_t width_right = matrix_right.width;

  size_t left_count = matrix_left.Count();
  size_t num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  std::vector<LocalMap> local_buffers(num_threads);
  std::vector<std::thread> threads(num_threads);

  size_t chunk = (left_count + num_threads - 1) / num_threads;

  for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t begin = thread_idx * chunk;
    size_t end = std::min(begin + chunk, left_count);

    threads[thread_idx] = std::thread([&, thread_idx, begin, end]() {
      ProcessChunk(begin, end, matrix_right, val_left, row_ind_left, col_ptr_left, local_buffers[thread_idx]);
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  std::map<Key, Complex> buffer;
  for (const auto &local : local_buffers) {
    for (const auto &[key, value] : local) {
      buffer[key] += value;
    }
  }

  CCSMatrix matrix_res;
  matrix_res.width = width_right;
  matrix_res.height = height_left;

  matrix_res.val.reserve(buffer.size());
  matrix_res.row_ind.reserve(buffer.size());
  matrix_res.col_ptr.reserve(buffer.size());

  for (const auto &[key, value] : buffer) {
    matrix_res.val.push_back(value);
    matrix_res.row_ind.push_back(key.first);
    matrix_res.col_ptr.push_back(key.second);
  }

  GetOutput() = matrix_res;
  return true;
}

bool PotashnikMMatrixMultComplexSTL::PostProcessingImpl() {
  return true;
}

}  // namespace potashnik_m_matrix_mult_complex
