#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "util/include/perf_test_util.hpp"
#include "votincev_d_radixmerge_sort/common/include/common.hpp"
#include "votincev_d_radixmerge_sort/seq/include/ops_seq.hpp"

namespace votincev_d_radixmerge_sort {

class VotincevDRadixMergeSortRunPerfTestsProcesses : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  InType GetTestInputData() final {
    return input_data;
  }

 protected:
  InType input_data;
  OutType expected_res;

  void SetUp() override {
    // Большой размер вектора для замера производительности
    size_t vect_sz = 1000000;
    input_data.assign(vect_sz, 0);
    for (size_t i = 0; i < vect_sz; i++) {
      // Генерация псевдослучайных чисел
      input_data[i] = static_cast<int32_t>((i * 1103515245 + 12345) % 2147483647);
    }

    expected_res = input_data;
    std::ranges::sort(expected_res);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // 1,2... процессы не владеют нужным результатом (актуально для MPI)
    if (output_data.size() != expected_res.size()) {
      return true;
    }
    // 0й процесс (или последовательный) должен иметь отсортированный массив
    return output_data == expected_res;
  }
};

namespace {
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, VotincevDRadixMergeSortSEQ>(PPC_SETTINGS_votincev_d_radixmerge_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = VotincevDRadixMergeSortRunPerfTestsProcesses::CustomPerfTestName;

TEST_P(VotincevDRadixMergeSortRunPerfTestsProcesses, RunPerfModes) {
  ExecuteTest(GetParam());
}

INSTANTIATE_TEST_SUITE_P(RunPerf, VotincevDRadixMergeSortRunPerfTestsProcesses, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace votincev_d_radixmerge_sort
