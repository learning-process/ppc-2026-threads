#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "goriacheva_k_mult_sparse_complex_matrix_ccs/common/include/common.hpp"
#include "goriacheva_k_mult_sparse_complex_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace goriacheva_k_mult_sparse_complex_matrix_ccs {

using Json = nlohmann::json;

class GoriachevaKMultSparseComplexMatrixCcsFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    const std::string path = ppc::util::GetAbsoluteTaskPath("goriacheva_k_mult_sparse_complex_matrix_ccs", "tests.txt");

    std::ifstream file(path);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open tests file: " + path);
    }

    Json tests_json;
    file >> tests_json;

    auto test_param = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    int test_index = std::get<0>(test_param);

    const Json &test = tests_json.at(test_index);

    input_data_ = {ParseMatrix(test["input"]["A"]), ParseMatrix(test["input"]["B"])};

    expected_ = ParseMatrix(test["result"]);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  bool CheckTestOutputData(OutType &output) final {
    return CompareMatrices(output, expected_);
  }

 private:
  static SparseMatrixCCS ParseMatrix(const Json &j) {
    SparseMatrixCCS m;
    m.rows = j.at("rows");
    m.cols = j.at("cols");

    for (const auto &val : j.at("values")) {
      m.values.emplace_back(val[0], val[1]);
    }

    m.row_ind = j.at("row_ind").get<std::vector<int>>();
    m.col_ptr = j.at("col_ptr").get<std::vector<int>>();
    return m;
  }

  static bool CompareMatrices(const SparseMatrixCCS &a, const SparseMatrixCCS &b) {
    if (a.rows != b.rows || a.cols != b.cols) {
      return false;
    }
    if (a.col_ptr != b.col_ptr) {
      return false;
    }
    if (a.row_ind != b.row_ind) {
      return false;
    }
    if (a.values.size() != b.values.size()) {
      return false;
    }

    const double eps = 1e-9;
    for (size_t i = 0; i < a.values.size(); i++) {
      if (std::abs(a.values[i].real() - b.values[i].real()) > eps ||
          std::abs(a.values[i].imag() - b.values[i].imag()) > eps) {
        return false;
      }
    }
    return true;
  }

  InType input_data_;
  SparseMatrixCCS expected_{};
};

namespace {

TEST_P(GoriachevaKMultSparseComplexMatrixCcsFuncTests, FromFile) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 14> kTestParams = {
    {TestType{0, "case0"}, TestType{1, "case1"}, TestType{2, "case2"}, TestType{3, "case3"}, TestType{4, "case4"},
     TestType{5, "case5"}, TestType{6, "case6"}, TestType{7, "case7"}, TestType{8, "case8"}, TestType{9, "case9"},
     TestType{10, "case10"}, TestType{11, "case11"}, TestType{12, "case12"}, TestType{13, "case13"}}};

const auto kTestTasksList = ppc::util::AddFuncTask<GoriachevaKMultSparseComplexMatrixCcsSEQ, InType>(
    kTestParams, PPC_SETTINGS_goriacheva_k_mult_sparse_complex_matrix_ccs);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName =
    GoriachevaKMultSparseComplexMatrixCcsFuncTests::PrintFuncTestName<GoriachevaKMultSparseComplexMatrixCcsFuncTests>;

INSTANTIATE_TEST_SUITE_P(FileMatrixTests, GoriachevaKMultSparseComplexMatrixCcsFuncTests, kGtestValues, kTestName);

}  // namespace
}  // namespace goriacheva_k_mult_sparse_complex_matrix_ccs
