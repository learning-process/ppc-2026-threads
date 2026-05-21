#include <gtest/gtest.h>

#include <cstddef>
#include <random>

#include "titaev_m_sortirovka_betchera/common/include/common.hpp"
#include "titaev_m_sortirovka_betchera/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevBatcherRadixPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  static constexpr size_t kSize = 200000;
  InType input;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10000.0, 10000.0);

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

TEST_P(TitaevBatcherRadixPerfTests, RunPerformanceModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TitaevSortirovkaBetcheraSEQ>(PPC_SETTINGS_titaev_m_sortirovka_betchera);

const auto kValues = ppc::util::TupleToGTestValues(kPerfTasks);
const auto kNameGen = TitaevBatcherRadixPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceSortingTests, TitaevBatcherRadixPerfTests, kValues, kNameGen);

}  // namespace
}  // namespace titaev_m_sortirovka_betchera
