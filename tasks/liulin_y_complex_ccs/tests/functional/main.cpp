#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "liulin_y_complex_ccs/common/include/common.hpp"
#include "liulin_y_complex_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace liulin_y_complex_ccs {

static CCSMatrix TripletToCcsTest(int rows_count, int cols_count,
                                  const std::vector<std::tuple<int, int, std::complex<double>>> &triplets) {
  CCSMatrix result{};
  result.count_rows = rows_count;
  result.count_cols = cols_count;
  result.col_index.assign(static_cast<size_t>(cols_count) + 1, 0);

  auto sorted_triplets{triplets};
  std::ranges::sort(sorted_triplets, [](const auto &lhs, const auto &rhs) {
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

  for (int col_idx{0}; col_idx < cols_count; ++col_idx) {
    result.col_index[static_cast<size_t>(col_idx) + 1] += result.col_index[static_cast<size_t>(col_idx)];
  }
  return result;
}

static void ReadMatrixFromFile(std::ifstream &file_stream, CCSMatrix &matrix) {
  int rows_val{0};
  int cols_val{0};
  int nnz_val{0};
  if (!(file_stream >> rows_val >> cols_val >> nnz_val)) {
    return;
  }
  std::vector<std::tuple<int, int, std::complex<double>>> triplets;
  for (int idx{0}; idx < nnz_val; ++idx) {
    int row_idx{0};
    int col_idx{0};
    double real_part{0.0};
    double imag_part{0.0};
    if (file_stream >> row_idx >> col_idx >> real_part >> imag_part) {
      triplets.emplace_back(row_idx, col_idx, std::complex<double>{real_part, imag_part});
    }
  }
  matrix = TripletToCcsTest(rows_val, cols_val, triplets);
}

static std::vector<std::complex<double>> ComputeDenseReference(const CCSMatrix &mat_a, const CCSMatrix &mat_b) {
  int rows_total{mat_a.count_rows};
  int cols_total{mat_b.count_cols};
  std::vector<std::complex<double>> dense(static_cast<size_t>(rows_total) * cols_total, {0.0, 0.0});

  for (int col_idx{0}; col_idx < cols_total; ++col_idx) {
    for (size_t idx_b{static_cast<size_t>(mat_b.col_index[col_idx])};
         idx_b < static_cast<size_t>(mat_b.col_index[col_idx + 1]); ++idx_b) {
      int mid_idx{mat_b.row_index[idx_b]};
      std::complex<double> val_b{mat_b.values[idx_b]};

      for (size_t idx_a{static_cast<size_t>(mat_a.col_index[mid_idx])};
           idx_a < static_cast<size_t>(mat_a.col_index[mid_idx + 1]); ++idx_a) {
        int row_idx{mat_a.row_index[idx_a]};
        dense[(static_cast<size_t>(row_idx) * cols_total) + col_idx] += mat_a.values[idx_a] * val_b;
      }
    }
  }
  return dense;
}

class LiulinYComplexCcsFuncTestsFromFile : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &params) {
    return std::get<1>(params);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string filename{std::get<1>(params)};
    std::string abs_path{ppc::util::GetAbsoluteTaskPath("liulin_y_complex_ccs", "seq/" + filename)};

    std::ifstream file_stream(abs_path + ".txt");
    if (!file_stream.is_open()) {
      throw std::runtime_error("Cannot open test file: " + abs_path + ".txt");
    }

    ReadMatrixFromFile(file_stream, input_data_.first);
    ReadMatrixFromFile(file_stream, input_data_.second);
    file_stream.close();

    int rows_res{input_data_.first.count_rows};
    int cols_res{input_data_.second.count_cols};
    auto dense_data{ComputeDenseReference(input_data_.first, input_data_.second)};

    std::vector<std::tuple<int, int, std::complex<double>>> res_triplets;
    for (int col_idx{0}; col_idx < cols_res; ++col_idx) {
      for (int row_idx{0}; row_idx < rows_res; ++row_idx) {
        auto value{dense_data[(static_cast<size_t>(row_idx) * cols_res) + col_idx]};
        if (std::abs(value) > 1e-15) {
          res_triplets.emplace_back(row_idx, col_idx, value);
        }
      }
    }
    exp_output_ = TripletToCcsTest(rows_res, cols_res, res_triplets);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count_rows != exp_output_.count_rows || output_data.count_cols != exp_output_.count_cols) {
      return false;
    }
    if (output_data.col_index != exp_output_.col_index || output_data.row_index != exp_output_.row_index) {
      return false;
    }
    for (size_t idx{0}; idx < output_data.values.size(); ++idx) {
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
