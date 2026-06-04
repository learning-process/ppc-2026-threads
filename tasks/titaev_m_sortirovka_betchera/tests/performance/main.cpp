#include <gtest/gtest.h>

#include <algorithm>
#include <tuple>

#include "titaev_m_sortirovka_betchera/all/include/ops_all.hpp"
#include "titaev_m_sortirovka_betchera/common/include/common.hpp"
#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"
#include "titaev_m_sortirovka_betchera/seq/include/ops_seq.hpp"
#include "titaev_m_sortirovka_betchera/stl/include/ops_stl.hpp"
#include "titaev_m_sortirovka_betchera/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevSortirovkaBetcheraPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    constexpr int kSize = 1 << 20;
    input_data_.resize(kSize);
    for (int i = 0; i < kSize; i++) {
      input_data_[i] = static_cast<double>(kSize - i) - 0.5;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(TitaevSortirovkaBetcheraPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kSeqPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TitaevSortirovkaBetcheraSEQ>(PPC_SETTINGS_titaev_m_sortirovka_betchera);
const auto kOmpPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TitaevSortirovkaBetcheraOMP>(PPC_SETTINGS_titaev_m_sortirovka_betchera);
const auto kTbbPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TitaevSortirovkaBetcheraTBB>(PPC_SETTINGS_titaev_m_sortirovka_betchera);
const auto kStlPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TitaevSortirovkaBetcheraSTL>(PPC_SETTINGS_titaev_m_sortirovka_betchera);
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TitaevSortirovkaBetcheraALL>(PPC_SETTINGS_titaev_m_sortirovka_betchera);

const auto kPerfTasks = std::tuple_cat(kSeqPerfTasks, kOmpPerfTasks, kTbbPerfTasks, kStlPerfTasks, kAllPerfTasks);
const auto kGtestValues = ppc::util::TupleToGTestValues(kPerfTasks);
const auto kPerfTestName = TitaevSortirovkaBetcheraPerfTests::CustomPerfTestName;
INSTANTIATE_TEST_SUITE_P(RunModeTests, TitaevSortirovkaBetcheraPerfTests, kGtestValues, kPerfTestName);
}  // namespace
}  // namespace titaev_m_sortirovka_betchera
