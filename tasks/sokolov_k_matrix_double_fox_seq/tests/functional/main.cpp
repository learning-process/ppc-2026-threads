#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "sokolov_k_matrix_double_fox_seq/common/include/common.hpp"
#include "sokolov_k_matrix_double_fox_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace sokolov_k_matrix_double_fox_seq {

namespace {

std::vector<double> MakeFilled(int n, double val) {
  return std::vector<double>(static_cast<std::size_t>(n * n), val);
}

std::vector<double> MakeIdentity(int n) {
  std::vector<double> m(static_cast<std::size_t>(n * n), 0.0);
  for (int i = 0; i < n; i++) {
    m[i * n + i] = 1.0;
  }
  return m;
}

}  // namespace

class SokolovKMatrixDoubleFoxFuncTestsSeq : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int n = std::get<1>(params);
    int b = std::get<2>(params);
    input_data_ = std::make_tuple(n, b, std::get<3>(params), std::get<4>(params));
    expected_ = std::get<5>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_.size()) {
      return false;
    }
    for (std::size_t i = 0; i < expected_.size(); i++) {
      if (std::abs(output_data[i] - expected_[i]) > 1e-9) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_;
};

namespace {

TEST_P(SokolovKMatrixDoubleFoxFuncTestsSeq, FoxMatmul) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {
    std::make_tuple("1x1_b1", 1, 1, std::vector<double>{5.0}, std::vector<double>{3.0}, std::vector<double>{15.0}),

    std::make_tuple("2x2_b1", 2, 1, std::vector<double>{1, 2, 3, 4}, std::vector<double>{5, 6, 7, 8},
                    std::vector<double>{19, 22, 43, 50}),

    std::make_tuple("2x2_b2", 2, 2, std::vector<double>{1, 2, 3, 4}, std::vector<double>{5, 6, 7, 8},
                    std::vector<double>{19, 22, 43, 50}),

    std::make_tuple("4x4_b1_ones", 4, 1, MakeFilled(4, 2.0), MakeFilled(4, 3.0), MakeFilled(4, 24.0)),

    std::make_tuple("4x4_b2_ones", 4, 2, MakeFilled(4, 2.0), MakeFilled(4, 3.0), MakeFilled(4, 24.0)),

    std::make_tuple("4x4_b4_ones", 4, 4, MakeFilled(4, 2.0), MakeFilled(4, 3.0), MakeFilled(4, 24.0)),

    std::make_tuple("4x4_b2_identity", 4, 2, MakeIdentity(4), MakeFilled(4, 7.0), MakeFilled(4, 7.0)),

    std::make_tuple("4x4_b2_zero", 4, 2, MakeFilled(4, 0.0), MakeFilled(4, 5.0), MakeFilled(4, 0.0)),

    std::make_tuple("6x6_b2", 6, 2, MakeFilled(6, 1.0), MakeFilled(6, 2.0), MakeFilled(6, 12.0)),

    std::make_tuple("6x6_b3", 6, 3, MakeFilled(6, 1.0), MakeFilled(6, 2.0), MakeFilled(6, 12.0)),

    std::make_tuple("3x3_b1", 3, 1, std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9},
                    std::vector<double>{9, 8, 7, 6, 5, 4, 3, 2, 1},
                    std::vector<double>{30, 24, 18, 84, 69, 54, 138, 114, 90}),

    std::make_tuple("3x3_b3", 3, 3, std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9},
                    std::vector<double>{9, 8, 7, 6, 5, 4, 3, 2, 1},
                    std::vector<double>{30, 24, 18, 84, 69, 54, 138, 114, 90}),
};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<SokolovKMatrixDoubleFoxSEQ, InType>(
    kTestParam, PPC_SETTINGS_sokolov_k_matrix_double_fox_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SokolovKMatrixDoubleFoxFuncTestsSeq::PrintFuncTestName<SokolovKMatrixDoubleFoxFuncTestsSeq>;

INSTANTIATE_TEST_SUITE_P(FoxMatmulTests, SokolovKMatrixDoubleFoxFuncTestsSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sokolov_k_matrix_double_fox_seq
