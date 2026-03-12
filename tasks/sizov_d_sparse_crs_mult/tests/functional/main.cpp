#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "sizov_d_sparse_crs_mult/common/include/common.hpp"
#include "sizov_d_sparse_crs_mult/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sizov_d_sparse_crs_mult {

namespace {

CRSMatrix DenseToCRS(const std::vector<std::vector<double>> &dense) {
  CRSMatrix m;
  m.rows = dense.size();
  m.cols = dense.empty() ? 0 : dense[0].size();
  m.row_ptr.resize(m.rows + 1, 0);

  for (std::size_t i = 0; i < m.rows; ++i) {
    for (std::size_t j = 0; j < m.cols; ++j) {
      if (std::abs(dense[i][j]) > 1e-12) {
        m.values.push_back(dense[i][j]);
        m.col_indices.push_back(j);
      }
    }
    m.row_ptr[i + 1] = m.values.size();
  }

  return m;
}

bool CompareCRS(const CRSMatrix &a, const CRSMatrix &b) {
  if (a.rows != b.rows || a.cols != b.cols) {
    return false;
  }
  if (a.row_ptr != b.row_ptr) {
    return false;
  }
  if (a.col_indices != b.col_indices) {
    return false;
  }
  if (a.values.size() != b.values.size()) {
    return false;
  }

  for (std::size_t i = 0; i < a.values.size(); ++i) {
    if (std::abs(a.values[i] - b.values[i]) > 1e-10) {
      return false;
    }
  }

  return true;
}

}  // namespace

class SizovDSparseCRSMultFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    return std::get<0>(param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    const auto &dense_a = std::get<1>(params);
    const auto &dense_b = std::get<2>(params);
    const auto &dense_c = std::get<3>(params);

    input_data_ = std::make_tuple(DenseToCRS(dense_a), DenseToCRS(dense_b));
    expected_ = DenseToCRS(dense_c);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return CompareCRS(output_data, expected_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CRSMatrix expected_;
};

namespace {

TEST_P(SizovDSparseCRSMultFuncTests, SpGemmCrsSeq) {
  ExecuteTest(GetParam());
}

using DenseMatrix = std::vector<std::vector<double>>;

const std::array<TestType, 6> kTestParams = {
    std::make_tuple("Identity_2x2", DenseMatrix{{1, 0}, {0, 1}}, DenseMatrix{{5, 6}, {7, 8}},
                    DenseMatrix{{5, 6}, {7, 8}}),

    std::make_tuple("ZeroMatrix", DenseMatrix{{0, 0}, {0, 0}}, DenseMatrix{{1, 2}, {3, 4}},
                    DenseMatrix{{0, 0}, {0, 0}}),

    std::make_tuple("Simple_2x2", DenseMatrix{{1, 2}, {3, 4}}, DenseMatrix{{5, 6}, {7, 8}},
                    DenseMatrix{{19, 22}, {43, 50}}),

    std::make_tuple("Rectangular_2x3_3x2", DenseMatrix{{1, 0, 2}, {0, 3, 0}}, DenseMatrix{{0, 1}, {4, 0}, {5, 6}},
                    DenseMatrix{{10, 13}, {12, 0}}),

    std::make_tuple("SparseDiagonal", DenseMatrix{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}},
                    DenseMatrix{{4, 0, 0}, {0, 5, 0}, {0, 0, 6}}, DenseMatrix{{4, 0, 0}, {0, 10, 0}, {0, 0, 18}}),

    std::make_tuple("SingleElement", DenseMatrix{{7}}, DenseMatrix{{8}}, DenseMatrix{{56}}),
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<SizovDSparseCRSMultSEQ, InType>(kTestParams, PPC_SETTINGS_sizov_d_sparse_crs_mult));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SizovDSparseCRSMultFuncTests::PrintFuncTestName<SizovDSparseCRSMultFuncTests>;

INSTANTIATE_TEST_SUITE_P(SparseCRSMultSEQTests, SizovDSparseCRSMultFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sizov_d_sparse_crs_mult
