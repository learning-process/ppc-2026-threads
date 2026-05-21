#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "makoveeva_matmul_double_stl/common/include/common.hpp"
#include "makoveeva_matmul_double_stl/stl/include/ops_stl.hpp"
#include "util/include/perf_test_util.hpp"

namespace makoveeva_matmul_double_stl {

class MatmulDoubleSTLPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  // Увеличиваем размер матрицы с 512 до 800 для более показательных результатов
  const int kMatrixSize_ = 800;
  InType input_data_;
  OutType expected_output_;

 protected:
  void SetUp() override {
    const auto n = static_cast<size_t>(kMatrixSize_);
    const size_t size = n * n;

    std::vector<double> a(size, 1.5);
    std::vector<double> b(size, 2.0);

    input_data_ = std::make_tuple(n, a, b);
    expected_output_.assign(size, 3.0 * static_cast<double>(n));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }

    const double epsilon = 1e-8;
    for (size_t i = 0; i < expected_output_.size(); ++i) {
      if (std::abs(expected_output_[i] - output_data[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(MatmulDoubleSTLPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, MatmulDoubleSTLTask>(PPC_SETTINGS_makoveeva_matmul_double_stl);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MatmulDoubleSTLPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, MatmulDoubleSTLPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace makoveeva_matmul_double_stl
