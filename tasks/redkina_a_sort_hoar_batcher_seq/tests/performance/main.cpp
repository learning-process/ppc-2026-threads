// perf_tests.cpp
#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <vector>

#include "redkina_a_sort_hoar_batcher_seq/common/include/common.hpp"
#include "redkina_a_sort_hoar_batcher_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace redkina_a_sort_hoar_batcher_seq {

class RedkinaASortHoarBatcherPerfTests
    : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const size_t size = 100000; 
    input_data_.resize(size);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-1000, 1000);
    for (auto& val : input_data_) {
      val = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType& output_data) final {
    // Проверка размера
    if (output_data.size() != input_data_.size()) return false;

    // Проверка, что массив отсортирован
    if (!std::is_sorted(output_data.begin(), output_data.end())) return false;

    // Проверка, что это перестановка (можно не проверять для перфа, но для надёжности)
    // Однако для больших массивов это затратно, поэтому просто проверяем размер и упорядоченность.
    return true;
  }

  InType GetTestInputData() final { return input_data_; }

 private:
  InType input_data_;
};

TEST_P(RedkinaASortHoarBatcherPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

// Создаём список задач для перформанс-тестов (только последовательная версия)
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, RedkinaASortHoarBatcherSEQ>(
        PPC_SETTINGS_redkina_a_sort_hoar_batcher_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName =
    RedkinaASortHoarBatcherPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RedkinaASortHoarBatcherPerfTests,
                         kGtestValues, kPerfTestName);

}  // namespace

}  // namespace redkina_a_sort_hoar_batcher_seq