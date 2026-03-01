#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <utility>
#include <vector>

#include "posternak_a_crs_mul_complex_matrix/common/include/common.hpp"
#include "posternak_a_crs_mul_complex_matrix/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace posternak_a_crs_mul_complex_matrix {

CRSMatrix MakeBandedCRS(int size, int bandwidth) {
  CRSMatrix m;
  m.rows = size;
  m.cols = size;
  m.index_row.reserve(size + 1);

  m.index_row.push_back(0);

  for (int row = 0; row < size; ++row) {
    int col_start = std::max(0, row - bandwidth);
    int col_end = std::min(size - 1, row + bandwidth);

    for (int col = col_start; col <= col_end; ++col) {
      auto real = static_cast<double>(row + 1);
      auto imag = static_cast<double>(col + 1);
      m.values.emplace_back(real, imag);
      m.index_col.push_back(col);
    }
    m.index_row.push_back(static_cast<int>(m.values.size()));
  }
  return m;
}

std::complex<double> ComputeExpectedValue(const CRSMatrix &a, const CRSMatrix &b, int row, int col) {
  std::complex<double> result(0.0, 0.0);

  for (int idx_a = a.index_row[row]; idx_a < a.index_row[row + 1]; ++idx_a) {
    int k = a.index_col[idx_a];

    for (int idx_b = b.index_row[k]; idx_b < b.index_row[k + 1]; ++idx_b) {
      if (b.index_col[idx_b] == col) {
        result += a.values[idx_a] * b.values[idx_b];
        break;
      }
    }
  }
  return result;
}

bool CheckSingleElement(const CRSMatrix &result, const CRSMatrix &a, const CRSMatrix &b, int row, int col) {
  std::complex<double> expected = ComputeExpectedValue(a, b, row, col);

  bool found = false;
  for (int idx = result.index_row[row]; idx < result.index_row[row + 1]; ++idx) {
    if (result.index_col[idx] == col) {
      found = true;
      if (std::abs(result.values[idx] - expected) > 1e-12) {
        return false;
      }
      break;
    }
  }

  if (std::abs(expected) > 1e-12 && !found) {
    return false;
  }
  if (std::abs(expected) <= 1e-12 && found) {
    return false;
  }

  return true;
}

bool CheckKeyElements(const CRSMatrix &result, const CRSMatrix &a, const CRSMatrix &b) {
  int n = result.rows;

  const std::vector<std::pair<int, int>> key_positions = {
      {0, 0}, {0, n - 1}, {n / 2, n / 2}, {n - 1, 0}, {n - 1, n - 1}};

  return std::ranges::all_of(key_positions, [&](const auto &pos) {
    const auto &[row, col] = pos;
    if (row >= n || col >= n) {
      return true;
    }
    return CheckSingleElement(result, a, b, row, col);
  });
}

class PosternakARunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int matrix_size = 10000;
    const int bandwidth = 40;

    CRSMatrix a = MakeBandedCRS(matrix_size, bandwidth);
    CRSMatrix b = MakeBandedCRS(matrix_size, bandwidth);

    input_data_ = {a, b};

    matrix_a_ = a;
    matrix_b_ = b;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return CheckKeyElements(output_data, matrix_a_, matrix_b_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CRSMatrix matrix_a_;
  CRSMatrix matrix_b_;
};

TEST_P(PosternakARunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, PosternakACRSMulComplexMatrixSEQ>(
    PPC_SETTINGS_posternak_a_crs_mul_complex_matrix);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PosternakARunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTests, PosternakARunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace posternak_a_crs_mul_complex_matrix
