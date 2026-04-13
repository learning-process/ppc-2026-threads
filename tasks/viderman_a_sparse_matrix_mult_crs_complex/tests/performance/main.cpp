#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <tuple>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/omp/include/ops_omp.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {
namespace {

bool CheckOutput(const CRSMatrix &out, int rows, int cols, std::size_t min_nnz,
                 std::optional<std::size_t> exact_nnz = std::nullopt) {
  if (!out.IsValid() || out.rows != rows || out.cols != cols) {
    return false;
  }
  if (out.NonZeros() < min_nnz) {
    return false;
  }
  if (exact_nnz.has_value() && out.NonZeros() != exact_nnz.value()) {
    return false;
  }
  return true;
}

CRSMatrix BuildBandMatrix(int n, const Complex &value, int bandwidth = 5) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    for (int offset = 0; offset < bandwidth; ++offset) {
      const int col = i + offset;
      if (col < n) {
        m.col_indices.push_back(col);
        m.values.push_back(value);
      }
    }
    m.row_ptr[i + 1] = static_cast<int>(m.col_indices.size());
  }
  return m;
}

CRSMatrix BuildDiagonalMatrix(int n, const Complex &value) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    m.col_indices.push_back(i);
    m.values.push_back(value);
    m.row_ptr[i + 1] = i + 1;
  }
  return m;
}

CRSMatrix BuildScatteredMatrix(int n, const Complex &value, int step = 7) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    const int col = (i * step) % n;
    m.col_indices.push_back(col);
    m.values.push_back(value);
    m.row_ptr[i + 1] = i + 1;
  }
  return m;
}

CRSMatrix BuildBlockDiagonalMatrix(int n, int block_size, const Complex &value) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    const int block_start = (i / block_size) * block_size;
    const int block_end = std::min(block_start + block_size, n);
    for (int col = block_start; col < block_end; ++col) {
      m.col_indices.push_back(col);
      m.values.push_back(value);
    }
    m.row_ptr[i + 1] = static_cast<int>(m.col_indices.size());
  }
  return m;
}

}  // namespace

class VidermanPerfBandMatrix : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const int n = 1000;
    input_data_ = std::make_tuple(BuildBandMatrix(n, Complex(1.0, 1.0)), BuildBandMatrix(n, Complex(1.0, 0.0)));
    rows_ = n;
    cols_ = n;
    min_nnz_ = 1U;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return CheckOutput(output_data, rows_, cols_, min_nnz_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  int rows_ = 0;
  int cols_ = 0;
  std::size_t min_nnz_ = 0U;
};

class VidermanPerfDiagonalMatrix : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const int n = 5000;
    input_data_ =
        std::make_tuple(BuildDiagonalMatrix(n, Complex(2.0, 1.0)), BuildDiagonalMatrix(n, Complex(1.0, -1.0)));
    rows_ = n;
    cols_ = n;
    exact_nnz_ = static_cast<std::size_t>(n);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (!exact_nnz_.has_value()) {
      return false;
    }
    return CheckOutput(output_data, rows_, cols_, *exact_nnz_, exact_nnz_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  int rows_ = 0;
  int cols_ = 0;
  std::optional<std::size_t> exact_nnz_;
};

class VidermanPerfWideBandMatrix : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const int n = 2000;
    input_data_ = std::make_tuple(BuildBandMatrix(n, Complex(1.0, 0.5), 20), BuildBandMatrix(n, Complex(0.5, 1.0), 20));
    rows_ = n;
    cols_ = n;
    min_nnz_ = 1U;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return CheckOutput(output_data, rows_, cols_, min_nnz_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  int rows_ = 0;
  int cols_ = 0;
  std::size_t min_nnz_ = 0U;
};

class VidermanPerfScatteredMatrix : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const int n = 3000;
    input_data_ =
        std::make_tuple(BuildScatteredMatrix(n, Complex(1.0, 1.0), 7), BuildScatteredMatrix(n, Complex(1.0, -1.0), 11));
    rows_ = n;
    cols_ = n;
    min_nnz_ = 1U;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return CheckOutput(output_data, rows_, cols_, min_nnz_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  int rows_ = 0;
  int cols_ = 0;
  std::size_t min_nnz_ = 0U;
};

class VidermanPerfBlockDiagonalMatrix : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const int n = 1000;
    const int block_size = 10;
    input_data_ = std::make_tuple(BuildBlockDiagonalMatrix(n, block_size, Complex(1.0, 1.0)),
                                  BuildBlockDiagonalMatrix(n, block_size, Complex(2.0, -1.0)));
    rows_ = n;
    cols_ = n;
    exact_nnz_ = static_cast<std::size_t>(n * block_size);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (!exact_nnz_.has_value()) {
      return false;
    }
    return CheckOutput(output_data, rows_, cols_, *exact_nnz_, exact_nnz_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  int rows_ = 0;
  int cols_ = 0;
  std::optional<std::size_t> exact_nnz_;
};

class VidermanPerfRectangularBandMatrix : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const int m = 500;
    const int k = 2000;
    const int n_out = 500;

    CRSMatrix a(m, k);
    for (int i = 0; i < m; ++i) {
      for (int offset = 0; offset < 5; ++offset) {
        const int col = (i * (k / m)) + offset;
        if (col < k) {
          a.col_indices.push_back(col);
          a.values.emplace_back(1.0, 0.5);
        }
      }
      a.row_ptr[i + 1] = static_cast<int>(a.col_indices.size());
    }

    CRSMatrix b(k, n_out);
    for (int i = 0; i < k; ++i) {
      const int col = (i * n_out) / k;
      if (col < n_out) {
        b.col_indices.push_back(col);
        b.values.emplace_back(0.5, 1.0);
      }
      b.row_ptr[i + 1] = static_cast<int>(b.col_indices.size());
    }

    input_data_ = std::make_tuple(a, b);
    rows_ = m;
    cols_ = n_out;
    min_nnz_ = 1U;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return CheckOutput(output_data, rows_, cols_, min_nnz_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  int rows_ = 0;
  int cols_ = 0;
  std::size_t min_nnz_ = 0U;
};

TEST_P(VidermanPerfBandMatrix, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(VidermanPerfDiagonalMatrix, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(VidermanPerfWideBandMatrix, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(VidermanPerfScatteredMatrix, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(VidermanPerfBlockDiagonalMatrix, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(VidermanPerfRectangularBandMatrix, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, VidermanASparseMatrixMultCRSComplexSEQ, VidermanASparseMatrixMultCRSComplexOMP>(
        PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ppc::util::BaseRunPerfTests<InType, OutType>::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(BandMatrixPerf, VidermanPerfBandMatrix, kGtestValues, kPerfTestName);
INSTANTIATE_TEST_SUITE_P(DiagonalMatrixPerf, VidermanPerfDiagonalMatrix, kGtestValues, kPerfTestName);
INSTANTIATE_TEST_SUITE_P(WideBandMatrixPerf, VidermanPerfWideBandMatrix, kGtestValues, kPerfTestName);
INSTANTIATE_TEST_SUITE_P(ScatteredMatrixPerf, VidermanPerfScatteredMatrix, kGtestValues, kPerfTestName);
INSTANTIATE_TEST_SUITE_P(BlockDiagonalMatrixPerf, VidermanPerfBlockDiagonalMatrix, kGtestValues, kPerfTestName);
INSTANTIATE_TEST_SUITE_P(RectangularBandMatrixPerf, VidermanPerfRectangularBandMatrix, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
