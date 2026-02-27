#include <gtest/gtest.h>

#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "potashnik_m_matrix_mult_complex/common/include/common.hpp"
#include "potashnik_m_matrix_mult_complex/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace potashnik_m_matrix_mult_complex {

class PotashnikMMatrixMultComplexPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kCount = 10000;
  InType input_data_;

  void SetUp() override {
    std::vector<std::vector<Complex>> matrix_left(kCount, std::vector<Complex>(kCount, Complex(0.0, 0.0)));
    std::vector<std::vector<Complex>> matrix_right(kCount, std::vector<Complex>(kCount, Complex(0.0, 0.0)));

    for (size_t i = 0; i < kCount; ++i) {
      matrix_left[i][(i * 11U) % kCount] = Complex(static_cast<double>(i) * 10.0, static_cast<double>(i) * 11.0);
      matrix_left[i][(i * 222U) % kCount] = Complex(static_cast<double>(i) * 9.0, static_cast<double>(i) * 12.0);
      matrix_right[i][(i * 55U) % kCount] = Complex(static_cast<double>(i) * 20.0, static_cast<double>(i) * 21.0);
      matrix_right[i][(i * 444U) % kCount] = Complex(static_cast<double>(i) * 19.0, static_cast<double>(i) * 22.0);
    }

    CCSMatrix matrix_left_ccs(matrix_left);
    CCSMatrix matrix_right_ccs(matrix_right);

    input_data_ = std::make_tuple(matrix_left_ccs, matrix_right_ccs);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const CCSMatrix &matrix_left = std::get<0>(input_data_);
    const CCSMatrix &matrix_right = std::get<1>(input_data_);

    std::vector<Complex> val_left = matrix_left.val;
    std::vector<size_t> row_ind_left = matrix_left.row_ind;
    std::vector<size_t> col_ptr_left = matrix_left.col_ptr;
    size_t height_left = matrix_left.height;
    // size_t width_left = matrix_left.width;

    std::vector<Complex> val_right = matrix_right.val;
    std::vector<size_t> row_ind_right = matrix_right.row_ind;
    std::vector<size_t> col_ptr_right = matrix_right.col_ptr;
    // size_t height_right = matrix_right.height;
    size_t width_right = matrix_right.width;

    std::map<std::pair<size_t, size_t>, Complex> buffer;

    for (size_t i = 0; i < matrix_left.Count(); i++) {
      size_t row_left = row_ind_left[i];
      size_t col_left = col_ptr_left[i];
      Complex left_val = val_left[i];

      for (size_t j = 0; j < matrix_right.Count(); j++) {
        size_t row_right = row_ind_right[j];
        size_t col_right = col_ptr_right[j];
        Complex right_val = val_right[j];

        if (col_left == row_right) {
          buffer[{row_left, col_right}] += left_val * right_val;
        }
      }
    }

    CCSMatrix matrix_res;

    matrix_res.width = width_right;
    matrix_res.height = height_left;
    for (const auto &[key, value] : buffer) {
      matrix_res.val.push_back(value);
      matrix_res.row_ind.push_back(key.first);
      matrix_res.col_ptr.push_back(key.second);
    }

    return output_data.Сompare(matrix_res);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(PotashnikMMatrixMultComplexPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PotashnikMMatrixMultComplexSEQ>(PPC_SETTINGS_potashnik_m_matrix_mult_complex);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PotashnikMMatrixMultComplexPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PotashnikMMatrixMultComplexPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace potashnik_m_matrix_mult_complex
