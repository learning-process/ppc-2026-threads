#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "zavyalov_a_complex_sparse_matr_mult/common/include/common.hpp"
#include "zavyalov_a_complex_sparse_matr_mult/seq/include/ops_seq.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

class ZavyalovAComplSparseMatrMultFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::to_string(std::get<1>(test_param)) + "_" +
           std::to_string(std::get<2>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    size_t rows_a = std::get<0>(params);
    size_t cols_a_rows_b = std::get<1>(params);
    size_t cols_b = std::get<2>(params);

    std::vector<std::vector<Complex>> matr_a(rows_a);
    for (size_t i = 0; i < rows_a; ++i) {
      matr_a[i].resize(cols_a_rows_b);
      for (size_t j = 0; j < cols_a_rows_b; ++j) {
        matr_a[i][j] = Complex(static_cast<double>((i * 42303U + 4242U + j) % 7433U),
                               static_cast<double>((i * 403U + 42U + j) % 733U));
      }
    }

    std::vector<std::vector<Complex>> matr_b(cols_a_rows_b);
    for (size_t i = 0; i < cols_a_rows_b; ++i) {
      matr_b[i].resize(cols_b);
      for (size_t j = 0; j < cols_b; ++j) {
        matr_b[i][j] = Complex(static_cast<double>((i * 42303U + 4242U + j) % 7433U),
                               static_cast<double>((i * 403U + 42U + j) % 733U));
      }
    }

    SparseMatrix matr1(matr_a);
    SparseMatrix matr2(matr_b);

    input_data_ = std::make_tuple(matr1, matr2);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const SparseMatrix &matr1 = std::get<0>(input_data_);
    const SparseMatrix &matr2 = std::get<1>(input_data_);

    size_t rows_a = matr1.height;
    size_t cols_a_rows_b = matr1.width;
    size_t cols_b = matr2.width;

    std::vector<std::vector<Complex>> matr_a(rows_a, std::vector<Complex>(cols_a_rows_b, Complex(0.0, 0.0)));
    std::vector<std::vector<Complex>> matr_b(cols_a_rows_b, std::vector<Complex>(cols_b, Complex(0.0, 0.0)));

    for (size_t idx = 0; idx < matr1.Count(); ++idx) {
      size_t row = matr1.row_ind[idx];
      size_t col = matr1.col_ind[idx];
      Complex val = matr1.val[idx];
      if (row < rows_a && col < cols_a_rows_b) {
        matr_a[row][col] = val;
      }
    }

    for (size_t idx = 0; idx < matr2.Count(); ++idx) {
      size_t row = matr2.row_ind[idx];
      size_t col = matr2.col_ind[idx];
      Complex val = matr2.val[idx];
      if (row < cols_a_rows_b && col < cols_b) {
        matr_b[row][col] = val;
      }
    }

    std::vector<std::vector<Complex>> matr_c(rows_a, std::vector<Complex>(cols_b, Complex(0.0, 0.0)));

    for (size_t i = 0; i < rows_a; ++i) {
      for (size_t j = 0; j < cols_b; ++j) {
        for (size_t k = 0; k < cols_a_rows_b; ++k) {
          matr_c[i][j] += (matr_a[i][k] * matr_b[k][j]);
        }
      }
    }

    SparseMatrix expected(matr_c);
    if (expected.Count() != output_data.Count()) {
      return false;
    }

    std::map<std::pair<size_t, size_t>, Complex> output_map;
    for (size_t idx = 0; idx < output_data.Count(); ++idx) {
      output_map[{output_data.row_ind[idx], output_data.col_ind[idx]}] = output_data.val[idx];
    }

    for (size_t idx = 0; idx < expected.Count(); ++idx) {
      auto key = std::make_pair(expected.row_ind[idx], expected.col_ind[idx]);
      auto it = output_map.find(key);
      if (it == output_map.end()) {
        return false;
      }
      if (!(expected.val[idx] == it->second)) {
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
};

namespace {

TEST_P(ZavyalovAComplSparseMatrMultFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {
    std::make_tuple(2, 3, 5), std::make_tuple(1, 1, 1), std::make_tuple(3, 3, 3), std::make_tuple(4, 3, 5),
    std::make_tuple(5, 4, 6), std::make_tuple(6, 5, 4), std::make_tuple(9, 3, 4), std::make_tuple(3, 9, 2),
    std::make_tuple(1, 5, 3), std::make_tuple(4, 7, 2)};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ZavyalovAComplSparseMatrMultSEQ, InType>(
    kTestParam, PPC_SETTINGS_zavyalov_a_complex_sparse_matr_mult));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    ZavyalovAComplSparseMatrMultFuncTests::PrintFuncTestName<ZavyalovAComplSparseMatrMultFuncTests>;

INSTANTIATE_TEST_SUITE_P(FuncTests, ZavyalovAComplSparseMatrMultFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zavyalov_a_compl_sparse_matr_mult
