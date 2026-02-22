#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
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

    size_t n = std::get<0>(params);
    size_t m = std::get<1>(params);
    size_t k = std::get<2>(params);

    std::vector<std::vector<Complex>> matr_a(n);
    for (size_t i = 0; i < n; i++) {
      matr_a[i].resize(m);
      for (size_t j = 0; j < m; j++) {
        matr_a[i][j] = Complex((i * 42303u + 4242u + j) % 7433u, (i * 403u + 42u + j) % 733u);
      }
    }

    std::vector<std::vector<Complex>> matr_b(m);
    for (size_t i = 0; i < m; i++) {
      matr_b[i].resize(k);
      for (size_t j = 0; j < k; j++) {
        matr_b[i][j] = Complex((i * 42303u + 4242u + j) % 7433u, (i * 403u + 42u + j) % 733u);
      }
    }

    Sparse_matrix matr1(matr_a);
    Sparse_matrix matr2(matr_b);

    input_data_ = std::make_tuple(matr1, matr2);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    Sparse_matrix &matr1 = std::get<0>(input_data_);
    Sparse_matrix &matr2 = std::get<1>(input_data_);
    size_t n = matr1.height;
    size_t m = matr1.width;
    size_t k = matr2.width;
    std::vector<std::vector<Complex>> matr_a(n);
    for (size_t i = 0; i < n; i++) {
      matr_a[i].resize(m);
    }
    std::vector<std::vector<Complex>> matr_b(m);
    for (size_t i = 0; i < m; i++) {
      matr_b[i].resize(k);
    }

    for (size_t i = 0; i < matr1.count(); i++) {
      size_t row = matr1.row_ind[i];
      size_t col = matr1.col_ind[i];
      Complex val = matr1.val[i];
      matr_a[row][col] = val;
    }

    for (size_t i = 0; i < matr2.count(); i++) {
      size_t row = matr2.row_ind[i];
      size_t col = matr2.col_ind[i];
      Complex val = matr2.val[i];
      matr_b[row][col] = val;
    }
    std::vector<std::vector<Complex>> matr_c(n, std::vector<Complex>(k, Complex(0, 0)));

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < k; j++) {
        for (size_t p = 0; p < m; p++) {
          matr_c[i][j] += (matr_a[i][p] * matr_b[p][j]);
        }
      }
    }

    Sparse_matrix res(matr_c);
    if (res.count() != output_data.count()) {
      return false;
    }

    for (size_t i = 0; i < res.count(); i++) {
      bool ok = true;
      ok &= (res.row_ind[i] == output_data.row_ind[i]);
      ok &= (res.col_ind[i] == output_data.col_ind[i]);
      ok &= (res.val[i] == output_data.val[i]);
      if (!ok) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
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
