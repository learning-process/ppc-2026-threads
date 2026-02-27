#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

// Заменяем пути example_threads на ваше название задачи
#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"
// В будущем здесь будут:
// #include "spichek_d_radix_sort_for_integers_with_simple_merging/omp/include/ops_omp.hpp"
// #include "spichek_d_radix_sort_for_integers_with_simple_merging/tbb/include/ops_tbb.hpp"

#include "util/include/perf_test_util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class RadixSortRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  // Для теста производительности берем большое количество элементов
  const int kCount_ = 1000000;
  InType input_data_{};

  void SetUp() override {
    input_data_.resize(kCount_);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-100000, 100000);

    for (int i = 0; i < kCount_; ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Для перформанс-теста достаточно убедиться, что массив отсортирован
    return std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(RadixSortRunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

// Если у вас еще нет OMP/TBB/STL версий, замените их на RadixSortSEQ
// или закомментируйте соответствующие аргументы в шаблоне MakeAllPerfTasks.
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType,
                                RadixSortSEQ,  // Заглушка для ALL
                                RadixSortSEQ,  // Заглушка для OMP
                                RadixSortSEQ,  // SEQ
                                RadixSortSEQ,  // Заглушка для STL
                                RadixSortSEQ   // Заглушка для TBB
                                >(PPC_SETTINGS_spichek_d_radix_sort_for_integers_with_simple_merging);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RadixSortRunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RadixSortRunPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
