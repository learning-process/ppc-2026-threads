#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "klimovich_v_crs_complex_mat_mul/common/include/common.hpp"
#include "klimovich_v_crs_complex_mat_mul/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace klimovich_v_crs_complex_mat_mul {
namespace {

constexpr double kCompareTol = 1e-9;

CrsMatrix FromDense(const std::vector<std::vector<Cplx>> &dense) {
  const int rows = static_cast<int>(dense.size());
  const int cols = rows == 0 ? 0 : static_cast<int>(dense.front().size());
  CrsMatrix out(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      const Cplx &v = dense[i][j];
      if (std::abs(v.real()) > 1e-15 || std::abs(v.imag()) > 1e-15) {
        out.col_indices.push_back(j);
        out.data.push_back(v);
      }
    }
    out.row_offsets[i + 1] = static_cast<int>(out.data.size());
  }
  return out;
}

std::vector<std::vector<Cplx>> ToDense(const CrsMatrix &m) {
  std::vector<std::vector<Cplx>> dense(m.n_rows, std::vector<Cplx>(m.n_cols, Cplx(0.0, 0.0)));
  for (int i = 0; i < m.n_rows; ++i) {
    for (int j = m.row_offsets[i]; j < m.row_offsets[i + 1]; ++j) {
      dense[i][m.col_indices[j]] = m.data[j];
    }
  }
  return dense;
}

std::vector<std::vector<Cplx>> GenerateDense(int rows, int cols, std::uint32_t seed_base) {
  std::vector<std::vector<Cplx>> dense(rows, std::vector<Cplx>(cols, Cplx(0.0, 0.0)));
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      const auto u = static_cast<std::uint32_t>((i * 31) + (j * 17) + seed_base);
      if (u % 4U == 0U) {
        const auto re = static_cast<double>(static_cast<int>(u % 17U) - 8);
        const auto im = static_cast<double>(static_cast<int>((u / 17U) % 13U) - 6);
        dense[i][j] = Cplx(re, im);
      }
    }
  }
  return dense;
}

std::vector<std::vector<Cplx>> DenseMultiply(const std::vector<std::vector<Cplx>> &lhs,
                                             const std::vector<std::vector<Cplx>> &rhs) {
  const auto rows = static_cast<int>(lhs.size());
  const int inner = rows == 0 ? 0 : static_cast<int>(lhs.front().size());
  const int cols = inner == 0 ? 0 : static_cast<int>(rhs.front().size());
  std::vector<std::vector<Cplx>> out(rows, std::vector<Cplx>(cols, Cplx(0.0, 0.0)));
  for (int i = 0; i < rows; ++i) {
    for (int k = 0; k < inner; ++k) {
      if (lhs[i][k] == Cplx(0.0, 0.0)) {
        continue;
      }
      for (int j = 0; j < cols; ++j) {
        out[i][j] += lhs[i][k] * rhs[k][j];
      }
    }
  }
  return out;
}

bool CrsApproxEqualDense(const CrsMatrix &actual, const std::vector<std::vector<Cplx>> &expected) {
  const auto dense_actual = ToDense(actual);
  if (dense_actual.size() != expected.size()) {
    return false;
  }
  for (std::size_t i = 0; i < expected.size(); ++i) {
    if (dense_actual[i].size() != expected[i].size()) {
      return false;
    }
    for (std::size_t j = 0; j < expected[i].size(); ++j) {
      const auto diff = dense_actual[i][j] - expected[i][j];
      if (std::abs(diff.real()) > kCompareTol || std::abs(diff.imag()) > kCompareTol) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace

class KlimovichVCrsComplexFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &p) {
    return std::to_string(std::get<0>(p)) + "x" + std::to_string(std::get<1>(p)) + "x" + std::to_string(std::get<2>(p));
  }

 protected:
  void SetUp() override {
    const auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int rows_a = std::get<0>(params);
    const int cols_a_rows_b = std::get<1>(params);
    const int cols_b = std::get<2>(params);

    auto dense_a = GenerateDense(rows_a, cols_a_rows_b, 11U);
    auto dense_b = GenerateDense(cols_a_rows_b, cols_b, 47U);

    expected_ = DenseMultiply(dense_a, dense_b);
    input_data_ = std::make_tuple(FromDense(dense_a), FromDense(dense_b));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    if (ppc::util::IsUnderMpirun()) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    if (rank != 0) {
      return true;
    }
    return CrsApproxEqualDense(output_data, expected_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  std::vector<std::vector<Cplx>> expected_;
};

namespace {

TEST_P(KlimovichVCrsComplexFuncTests, MultiplySparseCrsComplex) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParams = {
    std::make_tuple(1, 1, 1),  std::make_tuple(2, 3, 4),  std::make_tuple(4, 4, 4), std::make_tuple(5, 7, 3),
    std::make_tuple(8, 8, 8),  std::make_tuple(3, 9, 5),  std::make_tuple(6, 2, 7), std::make_tuple(10, 6, 10),
    std::make_tuple(11, 3, 4), std::make_tuple(7, 12, 9),
};

const auto kTaskList = std::tuple_cat(ppc::util::AddFuncTask<KlimovichVCrsComplexMatMulSeq, InType>(
    kTestParams, PPC_SETTINGS_klimovich_v_crs_complex_mat_mul));

const auto kGtestValues = ppc::util::ExpandToValues(kTaskList);
const auto kPerfTestName = KlimovichVCrsComplexFuncTests::PrintFuncTestName<KlimovichVCrsComplexFuncTests>;

INSTANTIATE_TEST_SUITE_P(FuncTests, KlimovichVCrsComplexFuncTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace klimovich_v_crs_complex_mat_mul
