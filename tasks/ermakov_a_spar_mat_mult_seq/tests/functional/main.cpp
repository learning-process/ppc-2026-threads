#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "ermakov_a_spar_mat_mult_seq/common/include/common.hpp"
#include "ermakov_a_spar_mat_mult_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace ermakov_a_spar_mat_mult_seq {

namespace {

using DenseMatrix = std::vector<std::vector<std::complex<double>>>;

DenseMatrix MakeDense(int rows, int cols) {
  return {static_cast<size_t>(rows), std::vector<std::complex<double>>(static_cast<size_t>(cols), {0.0, 0.0})};
}

void MultiplyRow(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &c, int i, int n, int p) {
  for (int k = 0; k < n; ++k) {
    const auto a_ik = a[i][k];
    if (a_ik == std::complex<double>(0.0, 0.0)) {
      continue;
    }

    for (int j = 0; j < p; ++j) {
      const auto b_kj = b[k][j];
      if (b_kj != std::complex<double>(0.0, 0.0)) {
        c[i][j] += a_ik * b_kj;
      }
    }
  }
}

DenseMatrix DenseMul(const DenseMatrix &a, const DenseMatrix &b) {
  const int m = static_cast<int>(a.size());
  const int n = (m != 0) ? static_cast<int>(a[0].size()) : 0;
  const int n_b = static_cast<int>(b.size());
  const int p = (n_b != 0) ? static_cast<int>(b[0].size()) : 0;

  DenseMatrix c = MakeDense(m, p);

  if (n != n_b) {
    return c;
  }

  for (int i = 0; i < m; ++i) {
    MultiplyRow(a, b, c, i, n, p);
  }

  return c;
}

double ResolveDensity(const std::string &desc) {
  if (desc == "VerySparse") {
    return 0.05;
  }
  if (desc == "MediumSparse") {
    return 0.2;
  }
  if (desc == "Dense") {
    return 0.7;
  }
  return 0.2;
}

void FillFixed(DenseMatrix &a, DenseMatrix &b) {
  a[0][0] = {1.0, 0.0};
  a[0][2] = {2.0, 0.0};
  a[2][1] = {3.0, 0.0};

  b[0][1] = {4.0, 0.0};
  b[1][2] = {5.0, 0.0};
  b[2][0] = {6.0, 0.0};
}

void FillRandom(DenseMatrix &a, DenseMatrix &b, int n, double density, std::mt19937 &gen) {
  std::uniform_real_distribution<double> dis_val(-5.0, 5.0);
  std::uniform_real_distribution<double> dis_prob(0.0, 1.0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (dis_prob(gen) < density) {
        a[i][j] = {dis_val(gen), dis_val(gen)};
      }
      if (dis_prob(gen) < density) {
        b[i][j] = {dis_val(gen), dis_val(gen)};
      }
    }
  }
}

MatrixCRS DenseToCRS(const DenseMatrix &m, double eps = 1e-12) {
  MatrixCRS r;

  const int rows = static_cast<int>(m.size());
  const int cols = (rows != 0) ? static_cast<int>(m[0].size()) : 0;

  r.rows = rows;
  r.cols = cols;
  r.row_ptr.resize(static_cast<size_t>(rows) + 1);
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

bool EqualCRS(const MatrixCRS &a, const MatrixCRS &b, double eps = 1e-9) {
  if (a.rows != b.rows || a.cols != b.cols) {
    return false;
  }
  if (a.row_ptr != b.row_ptr) {
    return false;
  }
  if (a.col_index != b.col_index) {
    return false;
  }
  if (a.values.size() != b.values.size()) {
    return false;
  }

  for (size_t i = 0; i < a.values.size(); ++i) {
    const auto &av = a.values[i];
    const auto &bv = b.values[i];

    if (std::abs(av.real() - bv.real()) > eps) {
      return false;
    }
    if (std::abs(av.imag() - bv.imag()) > eps) {
      return false;
    }
  }

  return true;
}

}  // namespace

class ErmakovARunFuncTestSparMatMult : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    const int n = std::get<0>(params);
    const std::string &desc = std::get<1>(params);

    DenseMatrix a = MakeDense(n, n);
    DenseMatrix b = MakeDense(n, n);

    if (desc == "SmallFixed") {
      if (n != 3) {
        GTEST_SKIP() << "SmallFixed defined only for n = 3";
      }
      FillFixed(a, b);
    } else {
      std::mt19937 gen(std::random_device{}());
      const double density = ResolveDensity(desc);
      FillRandom(a, b, n, density, gen);
    }

    input_data_.A = DenseToCRS(a);
    input_data_.B = DenseToCRS(b);

    DenseMatrix c_dense = DenseMul(a, b);
    expected_output_ = DenseToCRS(c_dense);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return EqualCRS(expected_output_, output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType expected_output_{};
};

namespace {

TEST_P(ErmakovARunFuncTestSparMatMult, MatmulCRSSeq) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {std::make_tuple(3, "SmallFixed"), std::make_tuple(10, "VerySparse"),
                                            std::make_tuple(20, "MediumSparse"), std::make_tuple(30, "Dense")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<ErmakovASparMatMultSEQ, InType>(kTestParam, PPC_SETTINGS_ermakov_a_spar_mat_mult_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ErmakovARunFuncTestSparMatMult::PrintFuncTestName<ErmakovARunFuncTestSparMatMult>;

INSTANTIATE_TEST_SUITE_P(SparseCRSMatMulTests, ErmakovARunFuncTestSparMatMult, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace ermakov_a_spar_mat_mult_seq
