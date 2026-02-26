#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include "liulin_y_complex_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace liulin_y_complex_ccs {

static CCSMatrix triplet_to_ccs_test(int rows, int cols,
                                     const std::vector<std::tuple<int, int, std::complex<double>>> &triplets) {
  CCSMatrix res;
  res.count_rows = rows;
  res.count_cols = cols;
  res.col_index.assign(cols + 1, 0);

  auto sorted = triplets;
  std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
    if (std::get<1>(a) != std::get<1>(b)) {
      return std::get<1>(a) < std::get<1>(b);
    }
    return std::get<0>(a) < std::get<0>(b);
  });

  for (const auto &triplet : sorted) {
    res.values.push_back(std::get<2>(triplet));
    res.row_index.push_back(std::get<0>(triplet));
    res.col_index[std::get<1>(triplet) + 1]++;
  }

  for (int i = 0; i < cols; ++i) {
    res.col_index[i + 1] += res.col_index[i];
  }
  return res;
}

class LiulinYComplexCcsFuncTestsFromFile : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &p) {
    return std::get<1>(p);
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

    auto read_matrix = [&](CCSMatrix &M) {
      int r;
      int c;
      int nnz;
      if (!(file >> r >> c >> nnz)) {
        return;
      }
      std::vector<std::tuple<int, int, std::complex<double>>> triplets;
      for (int i = 0; i < nnz; ++i) {
        int row;
        int col;
        double re;
        double im;
        if (file >> row >> col >> re >> im) {
          triplets.emplace_back(row, col, std::complex<double>(re, im));
        }
      }
      M = triplet_to_ccs_test(r, c, triplets);
    };

    read_matrix(input_data_.first);
    read_matrix(input_data_.second);

    int R = input_data_.first.count_rows;
    int C = input_data_.second.count_cols;
    std::vector<std::complex<double>> dense(R * C, {0.0, 0.0});

    for (int j = 0; j < C; ++j) {
      for (int kb = input_data_.second.col_index[j]; kb < input_data_.second.col_index[j + 1]; ++kb) {
        int k = input_data_.second.row_index[kb];
        std::complex<double> b_val = input_data_.second.values[kb];
        for (int ka = input_data_.first.col_index[k]; ka < input_data_.first.col_index[k + 1]; ++ka) {
          int i = input_data_.first.row_index[ka];
          dense[i * C + j] += input_data_.first.values[ka] * b_val;
        }
      }
    }

    std::vector<std::tuple<int, int, std::complex<double>>> res_triplets;
    for (int j = 0; j < C; ++j) {
      for (int i = 0; i < R; ++i) {
        if (std::abs(dense[i * C + j]) > 1e-15) {
          res_triplets.emplace_back(i, j, dense[i * C + j]);
        }
      }
    }
    exp_output_ = triplet_to_ccs_test(R, C, res_triplets);
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
    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i] - exp_output_.values[i]) > 1e-9) {
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
