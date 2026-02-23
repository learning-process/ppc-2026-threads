#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <random>
#include <vector>

#include "ermakov_a_spar_mat_mult_seq/common/include/common.hpp"
#include "ermakov_a_spar_mat_mult_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace ermakov_a_spar_mat_mult_seq {

namespace {

using DenseMatrix = std::vector<std::vector<std::complex<double>>>;

DenseMatrix MakeRandomDense(int n, double density) {
  DenseMatrix m(n, std::vector<std::complex<double>>(n, {0.0, 0.0}));

  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<double> dis_val(-5.0, 5.0);
  std::uniform_real_distribution<double> dis_prob(0.0, 1.0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (dis_prob(gen) < density) {
        m[i][j] = {dis_val(gen), dis_val(gen)};
      }
    }
  }
  return m;
}

MatrixCRS DenseToCRS(const DenseMatrix &m, double eps = 1e-12) {
  MatrixCRS r;

  const int rows = static_cast<int>(m.size());
  const int cols = (rows != 0) ? static_cast<int>(m[0].size()) : 0;

  r.rows = rows;
  r.cols = cols;
  r.row_ptr.resize(rows + 1);
  r.row_ptr[0] = 0;

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      const auto v = m[i][j];
      if (std::abs(v.real()) > eps || std::abs(v.imag()) > eps) {
        r.values.push_back(v);
        r.col_index.push_back(j);
      }
    }
    r.row_ptr[i + 1] = static_cast<int>(r.values.size());
  }

  return r;
}

}  // namespace

class ErmakovARunPerfTestSparMatMult : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int k_count = 15000;
  InType input_data{};

  void SetUp() override {
    const double density = 0.001;

    auto a_dense = MakeRandomDense(k_count, density);
    auto b_dense = MakeRandomDense(k_count, density);

    input_data.A = DenseToCRS(a_dense);
    input_data.B = DenseToCRS(b_dense);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != k_count || output_data.cols != k_count) {
      return false;
    }

    if (output_data.row_ptr.size() != static_cast<size_t>(output_data.rows) + 1) {
      return false;
    }

    if (output_data.values.size() != output_data.col_index.size()) {
      return false;
    }

    if (output_data.row_ptr.front() != 0) {
      return false;
    }

    if (static_cast<size_t>(output_data.row_ptr.back()) != output_data.values.size()) {
      return false;
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(ErmakovARunPerfTestSparMatMult, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ErmakovASparMatMultSEQ>(PPC_SETTINGS_ermakov_a_spar_mat_mult_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ErmakovARunPerfTestSparMatMult::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ErmakovARunPerfTestSparMatMult, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace ermakov_a_spar_mat_mult_seq
