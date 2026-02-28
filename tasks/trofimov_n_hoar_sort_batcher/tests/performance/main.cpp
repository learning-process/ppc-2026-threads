#include <cstddef>
#include <vector>
#include <random>
#include <algorithm>
#include <gtest/gtest.h>

#include "trofimov_n_hoar_sort_batcher/common/include/common.hpp"
#include "trofimov_n_hoar_sort_batcher/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace trofimov_n_hoar_sort_batcher {

using InType = std::vector<int>;
using OutType = std::vector<int>;

class TrofimovNHoarSortBatcherPerfTests :
    public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
    void SetUp() override {
        const size_t count = 300000;
        input_data_.resize(count);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dist(-1000000, 1000000);

        for(size_t i=0;i<count;++i) input_data_[i]=dist(gen);
    }

    bool CheckTestOutputData(OutType& output_data) final {
        return !output_data.empty() && std::ranges::is_sorted(output_data);
    }

    InType GetTestInputData() final { return input_data_; }

 private:
    InType input_data_;
};


TEST_P(TrofimovNHoarSortBatcherPerfTests, RunPerfTests) {
    ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, TrofimovNHoarSortBatcherSEQ>(
    PPC_SETTINGS_trofimov_n_hoar_sort_batcher);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = TrofimovNHoarSortBatcherPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(
    RunModeTests,
    TrofimovNHoarSortBatcherPerfTests,
    kGtestValues,
    kPerfTestName
);
} // namespace

} // namespace trofimov_n_hoar_sort_batcher