#include <gtest/gtest.h>

#include <cstddef>
#include <random>

#include "titaev_m_sortirovka_betchera/common/include/common.hpp"
#include "titaev_m_sortirovka_betchera/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevBatcherRadixPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  static constexpr size_t kSize = 1000000;
  InType input;

  void SetUp() override {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> dist(-100000.0, 100000.0);

    input.resize(kSize);
    for (auto &val : input) {
      val = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output) final {
    if (output.size() != input.size()) {
      return false;
    }
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

TEST_P(TitaevBatcherRadixPerfTests, RunPerformanceTBB) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TitaevSortirovkaBetcheraTBB>(PPC_SETTINGS_titaev_m_sortirovka_betchera);

const auto kValues = ppc::util::TupleToGTestValues(kPerfTasks);
const auto kNameGen = TitaevBatcherRadixPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceSortingTests, TitaevBatcherRadixPerfTests, kValues, kNameGen);

}  // namespace
}  // namespace titaev_m_sortirovka_betchera
