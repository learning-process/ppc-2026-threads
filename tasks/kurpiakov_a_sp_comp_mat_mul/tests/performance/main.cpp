#include <gtest/gtest.h>

#include <utility>
#include <vector>

#include "kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp"
#include "kurpiakov_a_sp_comp_mat_mul/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

namespace {

SparseMatrix MakeTridiagonal(int n) {
  SparseMatrix m;
  m._rows = n;
  m._cols = n;
  m._row_ptr.resize(n + 1, 0);

  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      m._values.emplace_back(0, 1);
      m._col_indices.push_back(i - 1);
    }
    m._values.emplace_back(2, 1);
    m._col_indices.push_back(i);
    if (i < n - 1) {
      m._values.emplace_back(1, 0);
      m._col_indices.push_back(i + 1);
    }
    m._row_ptr[i + 1] = static_cast<int>(m._values.size());
  }
  return m;
}

}  // namespace

class KurpiakovRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType expected_output_;

  void SetUp() override {
    constexpr int kSize = 500;

    auto a = MakeTridiagonal(kSize);
    auto b = MakeTridiagonal(kSize);

    expected_output_ = a.Multiply(b);
    input_data_ = {std::move(a), std::move(b)};
  }

  bool CheckTestOutputData(OutType& output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KurpiakovRunPerfTests, SparseMatMulPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KurpiskovACRSMatMulSEQ>(PPC_SETTINGS_kurpiakov_a_sp_comp_mat_mul);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = KurpiakovRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SparseMatMulPerfTests, KurpiakovRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kurpiakov_a_sp_comp_mat_mul
