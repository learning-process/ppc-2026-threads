#include "makoveeva_matmul_double_all/all/include/ops_all.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "makoveeva_matmul_double_all/common/include/common.hpp"

namespace makoveeva_matmul_double_all {
namespace {

constexpr size_t kSmallMatrixThreshold = 64;
constexpr size_t kMediumMatrixThreshold = 256;
constexpr size_t kLargeMatrixThreshold = 512;

void process_block_seq(const std::vector<double>& a,
                       const std::vector<double>& b,
                       std::vector<double>& c,
                       int n,
                       int i_start, int i_end,
                       int j_start, int j_end,
                       int k_start, int k_end) {
  for (int i = i_start; i < i_end; ++i) {
    for (int j = j_start; j < j_end; ++j) {
      double sum = 0.0;
      for (int k = k_start; k < k_end; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      c[(i * n) + j] += sum;
    }
  }
}

int calculate_block_size(int n) {
  return std::max(1, static_cast<int>(std::sqrt(static_cast<double>(n))));
}

int calculate_num_blocks(int n, int block_size) {
  return (n + block_size - 1) / block_size;
}

}  // namespace

MatmulDoubleAllTask::MatmulDoubleAllTask(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool MatmulDoubleAllTask::ValidationImpl() {
  const auto& input = GetInput();
  const size_t n = std::get<0>(input);
  const auto& a = std::get<1>(input);
  const auto& b = std::get<2>(input);

  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool MatmulDoubleAllTask::PreProcessingImpl() {
  const auto& input = GetInput();
  n_ = std::get<0>(input);
  a_ = std::get<1>(input);
  b_ = std::get<2>(input);
  c_.assign(n_ * n_, 0.0);

  choose_implementation();
  return true;
}

void MatmulDoubleAllTask::choose_implementation() {
  if (n_ <= kSmallMatrixThreshold) {
    selected_tech_ = TechType::kSeq;
  } else if (n_ <= kMediumMatrixThreshold) {
#ifdef _OPENMP
    selected_tech_ = TechType::kOmp;
#else
    selected_tech_ = TechType::kSeq;
#endif
  } else if (n_ <= kLargeMatrixThreshold) {
    selected_tech_ = TechType::kTbb;
  } else {
    selected_tech_ = TechType::kStl;
  }
}

void MatmulDoubleAllTask::multiply_seq() {
  if (n_ <= 0) return;

  c_.assign(c_.size(), 0.0);

  const int n_int = static_cast<int>(n_);
  const int block_size = calculate_block_size(n_int);
  const int num_blocks = calculate_num_blocks(n_int, block_size);

  for (int ib = 0; ib < num_blocks; ++ib) {
    for (int jb = 0; jb < num_blocks; ++jb) {
      for (int kb = 0; kb < num_blocks; ++kb) {
        const int i_start = ib * block_size;
        const int i_end = std::min(i_start + block_size, n_int);
        const int j_start = jb * block_size;
        const int j_end = std::min(j_start + block_size, n_int);
        const int k_start = kb * block_size;
        const int k_end = std::min(k_start + block_size, n_int);

        process_block_seq(a_, b_, c_, n_int, i_start, i_end, j_start, j_end, k_start, k_end);
      }
    }
  }
}

void MatmulDoubleAllTask::multiply_omp() {
#ifdef _OPENMP
  const size_t n = n_;
  const auto& a = a_;
  const auto& b = b_;
  auto& c = c_;

  #pragma omp parallel for collapse(2) default(none) shared(a, b, c, n)
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      c[(i * n) + j] = sum;
    }
  }
#else
  multiply_seq();
#endif
}

void MatmulDoubleAllTask::multiply_tbb() {
  const size_t n = n_;
  const auto& a = a_;
  const auto& b = b_;
  auto& c = c_;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
    [&](const tbb::blocked_range<size_t>& range) {
      for (size_t i = range.begin(); i < range.end(); ++i) {
        for (size_t j = 0; j < n; ++j) {
          double sum = 0.0;
          for (size_t k = 0; k < n; ++k) {
            sum += a[(i * n) + k] * b[(k * n) + j];
          }
          c[(i * n) + j] = sum;
        }
      }
    });
}

void MatmulDoubleAllTask::multiply_stl() {
  const size_t n = n_;
  const auto& a = a_;
  const auto& b = b_;
  auto& c = c_;

  const size_t num_threads = std::thread::hardware_concurrency();
  const size_t rows_per_thread = (n + num_threads - 1) / num_threads;

  std::vector<std::thread> threads;

  for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t start_row = thread_idx * rows_per_thread;
    size_t end_row = (std::min)(start_row + rows_per_thread, n);

    if (start_row >= n) break;

    threads.emplace_back([&a, &b, &c, n, start_row, end_row]() {
      for (size_t i = start_row; i < end_row; ++i) {
        for (size_t j = 0; j < n; ++j) {
          double sum = 0.0;
          for (size_t k = 0; k < n; ++k) {
            sum += a[(i * n) + k] * b[(k * n) + j];
          }
          c[(i * n) + j] = sum;
        }
      }
    });
  }

  for (auto& thread : threads) {
    thread.join();
  }
}

bool MatmulDoubleAllTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  switch (selected_tech_) {
    case TechType::kSeq:
      multiply_seq();
      break;
    case TechType::kOmp:
      multiply_omp();
      break;
    case TechType::kTbb:
      multiply_tbb();
      break;
    case TechType::kStl:
      multiply_stl();
      break;
    default:
      multiply_seq();
      break;
  }

  GetOutput() = c_;
  return true;
}

bool MatmulDoubleAllTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_all