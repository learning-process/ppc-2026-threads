#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../common/include/common.hpp"
#include "../../tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sabutay_sparse_complex_ccs_mult_tbb {

// Constants for test values
constexpr double kZero = 0.0;
constexpr double kValue1 = 1.0;
constexpr double kValue2 = 2.0;
constexpr double kValue3 = 3.0;
constexpr double kValue4 = 4.0;
constexpr double kValue5 = 5.0;
constexpr double kValue6 = 6.0;
constexpr double kValue12 = 12.0;
constexpr double kValue19 = 19.0;

class SabutayARunFuncTestsTbb : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  using TestTypeBase = typename std::remove_cv<typename std::remove_reference<TestType>::type>::type;
  static constexpr bool IsIntegralTestType = std::is_integral<TestTypeBase>::value;

  template <typename T = TestType>
  static auto PrintTestParam(const T &test_param) -> typename std::enable_if<
      std::is_integral<typename std::remove_cv<typename std::remove_reference<T>::type>::type>::value,
      std::string>::type {
    return std::to_string(test_param);
  }

  template <typename T = TestType>
  static auto PrintTestParam(const T &test_param) -> typename std::enable_if<
      !std::is_integral<typename std::remove_cv<typename std::remove_reference<T>::type>::type>::value,
      std::string>::type {
    return std::to_string(std::get<0>(test_param));
  }

  template <typename T = TestType>
  static auto GetTestParamIndex(const T &test_param) -> typename std::enable_if<
      std::is_integral<typename std::remove_cv<typename std::remove_reference<T>::type>::type>::value, int>::type {
    return static_cast<int>(test_param);
  }

  template <typename T = TestType>
  static auto GetTestParamIndex(const T &test_param) -> typename std::enable_if<
      !std::is_integral<typename std::remove_cv<typename std::remove_reference<T>::type>::type>::value, int>::type {
    return static_cast<int>(std::get<0>(test_param));
  }

  template <typename T = TestType>
  static auto MakeTestParam(int idx) -> typename std::enable_if<
      std::is_integral<typename std::remove_cv<typename std::remove_reference<T>::type>::type>::value, T>::type {
    return static_cast<T>(idx);
  }

  template <typename T = TestType>
  static auto MakeTestParam(int idx) -> typename std::enable_if<
      !std::is_integral<typename std::remove_cv<typename std::remove_reference<T>::type>::type>::value, T>::type {
    return T{idx, std::string{}};
  }

  template <typename Matrix>
  static auto MatrixRowsImpl(Matrix &matrix, int) -> decltype(matrix.m) & {
    return matrix.m;
  }

  template <typename Matrix>
  static auto MatrixRowsImpl(Matrix &matrix, long) -> int & {
    return matrix.rows;
  }

  template <typename Matrix>
  static auto MatrixRows(Matrix &matrix) -> int & {
    return MatrixRowsImpl(matrix, 0);
  }

  template <typename Matrix>
  static auto MatrixColsImpl(Matrix &matrix, int) -> decltype(matrix.n) & {
    return matrix.n;
  }

  template <typename Matrix>
  static auto MatrixColsImpl(Matrix &matrix, long) -> int & {
    return matrix.cols;
  }

  template <typename Matrix>
  static auto MatrixCols(Matrix &matrix) -> int & {
    return MatrixColsImpl(matrix, 0);
  }

  template <typename Matrix>
  static auto MatrixRowsImpl(const Matrix &matrix, int) -> decltype(matrix.m) {
    return matrix.m;
  }

  template <typename Matrix>
  static auto MatrixRowsImpl(const Matrix &matrix, long) -> int {
    return matrix.rows;
  }

  template <typename Matrix>
  static auto MatrixRows(const Matrix &matrix) -> int {
    return MatrixRowsImpl(matrix, 0);
  }

  template <typename Matrix>
  static auto MatrixColsImpl(const Matrix &matrix, int) -> decltype(matrix.n) {
    return matrix.n;
  }

  template <typename Matrix>
  static auto MatrixColsImpl(const Matrix &matrix, long) -> int {
    return matrix.cols;
  }

  template <typename Matrix>
  static auto MatrixCols(const Matrix &matrix) -> int {
    return MatrixColsImpl(matrix, 0);
  }

 protected:
  void SetUp() override {
    const auto params =
        GetTestParamIndex(std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));
    auto &a = std::get<0>(input_data_);
    auto &b = std::get<1>(input_data_);
    auto &c = test_result_;

    if (params == 0) {
      MatrixRows(a) = 2;
      MatrixCols(a) = 3;
      a.col_ptr = {0, 1, 2, 3};
      a.row_ind = {0, 1, 0};
      a.values = {{kValue1, kZero}, {kValue2, kZero}, {kValue3, kZero}};

      MatrixRows(b) = 3;
      MatrixCols(b) = 2;
      b.col_ptr = {0, 2, 3};
      b.row_ind = {0, 2, 1};
      b.values = {{kValue4, kZero}, {kValue5, kZero}, {kValue6, kZero}};

      MatrixRows(c) = 2;
      MatrixCols(c) = 2;
      c.col_ptr = {0, 1, 2};
      c.row_ind = {0, 1};
      c.values = {{kValue19, kZero}, {kValue12, kZero}};
    }
    if (params == 1) {
      MatrixRows(a) = 2;
      MatrixCols(a) = 2;
      a.col_ptr = {0, 1, 2};
      a.row_ind = {0, 1};
      a.values = {{kValue1, kZero}, {kValue2, kZero}};

      MatrixRows(b) = 2;
      MatrixCols(b) = 2;
      b.col_ptr = {0, 1, 2};
      b.row_ind = {1, 0};
      b.values = {{kValue3, kZero}, {kValue4, kZero}};

      MatrixRows(c) = 2;
      MatrixCols(c) = 2;
      c.col_ptr = {0, 1, 2};
      c.row_ind = {1, 0};
      c.values = {{kValue6, kZero}, {kValue4, kZero}};
    }
    if (params == 2) {
      MatrixRows(a) = 1;
      MatrixCols(a) = 1;
      a.col_ptr = {0, 1};
      a.row_ind = {0};
      a.values = {{kValue2, kValue1}};

      MatrixRows(b) = 1;
      MatrixCols(b) = 1;
      b.col_ptr = {0, 1};
      b.row_ind = {0};
      b.values = {{kValue1, kValue1}};

      MatrixRows(c) = 1;
      MatrixCols(c) = 1;
      c.col_ptr = {0, 1};
      c.row_ind = {0};
      c.values = {{kValue1, kValue3}};
    }
  }

  bool CheckTestOutputData(OutType &output_data) override {
    bool result = true;
    constexpr double kEps = 1e-14;
    if (MatrixRows(test_result_) != MatrixRows(output_data) || MatrixCols(test_result_) != MatrixCols(output_data) ||
        test_result_.col_ptr.size() != output_data.col_ptr.size() ||
        test_result_.row_ind.size() != output_data.row_ind.size() ||
        test_result_.values.size() != output_data.values.size()) {
      return false;
    }

    for (size_t i = 0; i < test_result_.col_ptr.size(); ++i) {
      if (test_result_.col_ptr[i] != output_data.col_ptr[i]) {
        return false;
      }
    }

    for (int j = 0; j < MatrixCols(test_result_); ++j) {
      std::vector<std::pair<int, std::complex<double>>> test;
      std::vector<std::pair<int, std::complex<double>>> output;
      for (int k = 0; k < test_result_.col_ptr[j + 1] - test_result_.col_ptr[j]; ++k) {
        test.emplace_back(test_result_.row_ind[test_result_.col_ptr[j] + k],
                          test_result_.values[test_result_.col_ptr[j] + k]);
        output.emplace_back(output_data.row_ind[output_data.col_ptr[j] + k],
                            output_data.values[output_data.col_ptr[j] + k]);
      }
      auto cmp = [](const auto &x, const auto &y) { return x.first < y.first; };
      std::sort(test.begin(), test.end(), cmp);
      std::sort(output.begin(), output.end(), cmp);
      for (size_t i = 0; i < test.size(); ++i) {
        if (test[i].first != output[i].first || std::abs(test[i].second - output[i].second) > kEps) {
          result = false;
        }
      }
    }

    return result;
  }

  InType GetTestInputData() override {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType test_result_;
};

namespace {

TEST_P(SabutayARunFuncTestsTbb, FuncCCSTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {SabutayARunFuncTestsTbb::MakeTestParam(0),
                                            SabutayARunFuncTestsTbb::MakeTestParam(1),
                                            SabutayARunFuncTestsTbb::MakeTestParam(2)};

const auto kTestTasksList = ppc::util::AddFuncTask<SabutayASparseComplexCcsMultTBB, InType>(
    kTestParam, PPC_SETTINGS_sabutay_sparse_complex_ccs_mult_tbb);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SabutayARunFuncTestsTbb::PrintFuncTestName<SabutayARunFuncTestsTbb>;

INSTANTIATE_TEST_SUITE_P(RunFuncCCSTest, SabutayARunFuncTestsTbb, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
