#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "titaev_m_sortirovka_betchera/common/include/common.hpp"
#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"
#include "titaev_m_sortirovka_betchera/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevBatcherRadixPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  static constexpr size_t kSize = 100000;
  InType input;
  void SetUp() override {
    std::mt19937 gen(123);
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
    input.resize(kSize);
    for (size_t i = 0; i < kSize; i++) {
      input[i] = dist(gen);
    }
  }
  bool CheckTestOutputData(OutType &output) final {
    for (size_t i = 1; i < output.size(); i++) {
      if (output[i] < output[i - 1]) {
        return false;
      }
    }
    return true;
  }
  InType GetTestInputData() final {
    return input;
  }
};

TEST_P(TitaevBatcherRadixPerfTests, RunPerf) {
  ExecuteTest(GetParam());
}

namespace {
inline std::string SafeGetSettings() {
  std::string s = PPC_SETTINGS_titaev_m_sortirovka_betchera;
  return (s.empty() || s == "null") ? "{}" : s;
}
const auto kPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TitaevSortirovkaBetcheraSEQ, TitaevSortirovkaBetcheraOMP>(SafeGetSettings());
INSTANTIATE_TEST_SUITE_P(PerformanceTests, TitaevBatcherRadixPerfTests, ppc::util::TupleToGTestValues(kPerfTasks),
                         TitaevBatcherRadixPerfTests::CustomPerfTestName);
}  // namespace
}  // namespace titaev_m_sortirovka_betchera
