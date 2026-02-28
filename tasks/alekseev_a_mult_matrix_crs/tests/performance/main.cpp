#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "alekseev_a_mult_matrix_crs/common/include/common.hpp"
#include "alekseev_a_mult_matrix_crs/seq/include/ops_seq.hpp"
#include "performance/include/performance.hpp"
#include "util/include/perf_test_util.hpp"

namespace alekseev_a_mult_matrix_crs {

namespace {

CRSMatrix GenerateBandMatrix(std::size_t n, std::size_t bandwidth, double value) {
  CRSMatrix m;
  m.rows = n;
  m.cols = n;
  m.row_ptr.resize(n + 1, 0);

  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t begin = (i > bandwidth ? i - bandwidth : 0);
    const std::size_t end = std::min(n - 1, i + bandwidth);
    for (std::size_t j = begin; j <= end; ++j) {
      m.values.push_back(value);
      m.col_indices.push_back(j);
    }
    m.row_ptr[i + 1] = m.values.size();
  }

  return m;
}

}  // namespace

class AlekseevAMultMatrixCRSPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  static constexpr std::size_t kSize = 40000;
  static constexpr std::size_t kBandwidth = 30;

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    ppc::util::BaseRunPerfTests<InType, OutType>::SetPerfAttributes(perf_attrs);
    // Measure one heavy run; this keeps wall time reasonable and metric readable.
    perf_attrs.num_running = 1;
  }

  void SetUp() override {
    CRSMatrix a = GenerateBandMatrix(kSize, kBandwidth, 2.0);
    CRSMatrix b = GenerateBandMatrix(kSize, kBandwidth, 3.0);
    input_data_ = std::make_tuple(a, b);
    // The exact product values are not used in perf test verification.
    expected_.rows = kSize;
    expected_.cols = kSize;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.rows == expected_.rows && output_data.cols == expected_.cols && !output_data.values.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CRSMatrix expected_;
};

TEST_P(AlekseevAMultMatrixCRSPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, AlekseevAMultMatrixCRSSEQ>(PPC_SETTINGS_alekseev_a_mult_matrix_crs);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = AlekseevAMultMatrixCRSPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(AlekseevAMultMatrixCRSPerfTests, AlekseevAMultMatrixCRSPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace alekseev_a_mult_matrix_crs
