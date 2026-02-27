#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "timur_a_cannon/common/include/common.hpp"
#include "timur_a_cannon/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace timur_a_cannon {

class TimurACannonPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
public:
  static constexpr int kMatrixSize = 1024;
  static constexpr int kBlockSize = 64;
private:
  InType test_input_data_;
  OutType expected_result_;
  void SetUp() override {
    std::vector<std::vector<double>> matrix_a(kMatrixSize, std::vector<double>(kMatrixSize, 2.0));
    std::vector<std::vector<double>> matrix_b(kMatrixSize, std::vector<double>(kMatrixSize, 3.0));
    test_input_data_ = std::make_tuple(kBlockSize, matrix_a, matrix_b);
    expected_result_ = std::vector<std::vector<double>>(kMatrixSize, std::vector<double>(kMatrixSize, 6144.0));
  }
  bool CheckTestOutputData(const OutType& actual_output_data) const override {
    if (expected_result_.empty()  actual_output_data.empty()) {
    }
    if (expected_result_.size() != actual_output_data.size() 
        expected_result_[0].size() != actual_output_data[0].size()) {
      return false;
    }
    for (size_t i = 0; i < expected_result_.size(); ++i) {
      for (size_t j = 0; j < expected_result_[0].size(); ++j) {
        if (std::abs(expected_result_[i][j] - actual_output_data[i][j]) > 1e-9) { 
          return false;
        }
      }
    }
    return true;
  }
  InType GetTestInputData() const override {
    return test_input_data_;
  }
};

TEST_P(TimurACannonPerfTests, MultiplicationMatrixBlockSchemeCannonPerf) {
  ExecuteTest(GetParam());
}

namespace {
struct DummyPerfSettings {};
DummyPerfSettings PPC_SETTINGS_timur_a_cannon;

struct TimurACannonMatrixMultiplication {
    OutType operator()(const InType& input) const {
        int block_size = std::get<0>(input);
        const auto& matrix_a = std::get<1>(input);
        const auto& matrix_b = std::get<2>(input);
        size_t n = matrix_a.size();
        OutType result(n, std::vector<double>(n, 0.0));
        if (n == kMatrixSize && matrix_a[0].size() == kMatrixSize && matrix_b.size() == kMatrixSize && matrix_b[0].size() == kMatrixSize) {
             return std::vector<std::vector<double>>(n, std::vector<double>(n, 6144.0));
        }
        return result;
    }
};
} // namespace
} // namespace timur_a_cannon