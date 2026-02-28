#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "alekseev_a_mult_matrix_crs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace alekseev_a_mult_matrix_crs {
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
  if (a.rows != b.rows || a.cols != b.cols || a.row_ptr != b.row_ptr || a.col_indices != b.col_indices) {
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

class AlekseevAMultMatrixCRSFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::make_tuple(DenseToCRS(std::get<1>(params)), DenseToCRS(std::get<2>(params)));
    expected_ = DenseToCRS(std::get<3>(params));
  }
  bool CheckTestOutputData(OutType &out) final {
    return CompareCRS(out, expected_);
  }
  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType expected_{};
};

TEST_P(AlekseevAMultMatrixCRSFuncTests, SeqRun) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kParams = {
    std::make_tuple("Identity", std::vector<std::vector<double>>{{1, 0}, {0, 1}},
                    std::vector<std::vector<double>>{{5, 6}, {7, 8}}, std::vector<std::vector<double>>{{5, 6}, {7, 8}}),
    std::make_tuple("Zero", std::vector<std::vector<double>>{{0, 0}, {0, 0}},
                    std::vector<std::vector<double>>{{1, 2}, {3, 4}}, std::vector<std::vector<double>>{{0, 0}, {0, 0}}),
    std::make_tuple("Simple", std::vector<std::vector<double>>{{1, 2}, {3, 4}},
                    std::vector<std::vector<double>>{{5, 6}, {7, 8}},
                    std::vector<std::vector<double>>{{19, 22}, {43, 50}}),
    std::make_tuple("Rect", std::vector<std::vector<double>>{{1, 0, 2}, {0, 3, 0}},
                    std::vector<std::vector<double>>{{0, 1}, {4, 0}, {5, 6}},
                    std::vector<std::vector<double>>{{10, 13}, {12, 0}})};

INSTANTIATE_TEST_SUITE_P(AlekseevAMultSEQ, AlekseevAMultMatrixCRSFuncTests,
                         ppc::util::ExpandToValues(ppc::util::AddFuncTask<AlekseevAMultMatrixCRSSEQ, InType>(
                             kParams, "alekseev_a_mult_matrix_crs")));

}  // namespace alekseev_a_mult_matrix_crs
