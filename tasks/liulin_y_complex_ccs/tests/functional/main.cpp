#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include "liulin_y_complex_ccs/common/include/common.hpp"
#include "liulin_y_complex_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace liulin_y_complex_ccs {

static CCSMatrix triplet_to_ccs_test(int rows, int cols,
                                     const std::vector<std::tuple<int, int, std::complex<double>>> &triplets) {
  CCSMatrix result;
  result.count_rows = rows;
  result.count_cols = cols;
  result.col_index.assign(static_cast<size_t>(cols) + 1, 0);

  auto sorted_triplets = triplets;
  std::sort(sorted_triplets.begin(), sorted_triplets.end(), [](const auto &lhs, const auto &rhs) {
    if (std::get<1>(lhs) != std::get<1>(rhs)) {
      return std::get<1>(lhs) < std::get<1>(rhs);
    }
    return std::get<0>(lhs) < std::get<0>(rhs);
  });

  for (const auto &triplet : sorted_triplets) {
    result.values.push_back(std::get<2>(triplet));
    result.row_index.push_back(std::get<0>(triplet));
    result.col_index[static_cast<size_t>(std::get<1>(triplet)) + 1]++;
  }

  for (int col_idx = 0; col_idx < cols; ++col_idx) {
    result.col_index[static_cast<size_t>(col_idx) + 1] += result.col_index[static_cast<size_t>(col_idx)];
  }
  return result;
}

class LiulinYComplexCcsFuncTestsFromFile : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &params) {
    return std::get<1>(params);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string filename = std::get<1>(params);
    std::string abs_path = ppc::util::GetAbsoluteTaskPath("liulin_y_complex_ccs", "seq/" + filename);

    std::ifstream file(abs_path + ".txt");
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open test file: " + abs_path + ".txt");
    }

    auto read_matrix = [&](CCSMatrix &matrix) {
      int rows_val = 0;
      int cols_val = 0;
      int nnz_val = 0;
      if (!(file >> rows_val >> cols_val >> nnz_val)) {
        return;
      }
      std::vector<std::tuple<int, int, std::complex<double>>> triplets;
      for (int idx = 0; idx < nnz_val; ++idx) {
        int row_idx = 0;
        int col_idx = 0;
        double real_part = 0.0;
        double imag_part = 0.0;
        if (file >> row_idx >> col_idx >> real_part >> imag_part) {
          triplets.emplace_back(row_idx, col_idx, std::complex<double>(real_part, imag_part));
        }
      }
      matrix = triplet_to_ccs_test(rows_val, cols_val, triplets);
    };

    read_matrix(input_data_.first);
    read_matrix(input_data_.second);

    int rows_total = input_data_.first.count_rows;
    int cols_total = input_data_.second.count_cols;
    std::vector<std::complex<double>> dense(static_cast<size_t>(rows_total) * cols_total, {0.0, 0.0});

    for (int col_idx = 0; col_idx < cols_total; ++col_idx) {
      for (int idx_b = input_data_.second.col_index[static_cast<size_t>(col_idx)];
           idx_b < input_data_.second.col_index[static_cast<size_t>(col_idx) + 1]; ++idx_b) {
        int mid_idx = input_data_.second.row_index[static_cast<size_t>(idx_b)];
        std::complex<double> val_b = input_data_.second.values[static_cast<size_t>(idx_b)];
        for (int idx_a = input_data_.first.col_index[static_cast<size_t>(mid_idx)];
             idx_a < input_data_.first.col_index[static_cast<size_t>(mid_idx) + 1]; ++idx_a) {
          int row_idx = input_data_.first.row_index[static_cast<size_t>(idx_a)];
          dense[(static_cast<size_t>(row_idx) * cols_total) + col_idx] +=
              input_data_.first.values[static_cast<size_t>(idx_a)] * val_b;
        }
      }
    }

    std::vector<std::tuple<int, int, std::complex<double>>> res_triplets;
    for (int col_idx = 0; col_idx < cols_total; ++col_idx) {
      for (int row_idx = 0; row_idx < rows_total; ++row_idx) {
        if (std::abs(dense[(static_cast<size_t>(row_idx) * cols_total) + col_idx]) > 1e-15) {
          res_triplets.emplace_back(row_idx, col_idx, dense[(static_cast<size_t>(row_idx) * cols_total) + col_idx]);
        }
      }
    }
    exp_output_ = triplet_to_ccs_test(rows_total, cols_total, res_triplets);
    file.close();
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count_rows != exp_output_.count_rows) {
      return false;
    }
    if (output_data.count_cols != exp_output_.count_cols) {
      return false;
    }
    if (output_data.col_index != exp_output_.col_index) {
      return false;
    }
    if (output_data.row_index != exp_output_.row_index) {
      return false;
    }
    if (output_data.values.size() != exp_output_.values.size()) {
      return false;
    }
    for (size_t idx = 0; idx < output_data.values.size(); ++idx) {
      if (std::abs(output_data.values[idx] - exp_output_.values[idx]) > 1e-9) {
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
  OutType exp_output_;
};

namespace {
TEST_P(LiulinYComplexCcsFuncTestsFromFile, SparseMultiplyFileTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParam = {
    std::make_tuple(0, "identity_2x2"),        std::make_tuple(1, "complex_scalar"),
    std::make_tuple(2, "rectangular_simple"),  std::make_tuple(3, "zero_matrix"),
    std::make_tuple(4, "sparse_random_small"), std::make_tuple(5, "only_imaginary")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<LiulinYComplexCcs, InType>(kTestParam, PPC_SETTINGS_liulin_y_complex_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kFuncTestName = LiulinYComplexCcsFuncTestsFromFile::PrintFuncTestName<LiulinYComplexCcsFuncTestsFromFile>;

INSTANTIATE_TEST_SUITE_P(Sequential, LiulinYComplexCcsFuncTestsFromFile, kGtestValues, kFuncTestName);
}  // namespace
}  // namespace liulin_y_complex_ccs
