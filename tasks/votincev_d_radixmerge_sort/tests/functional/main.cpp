#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "votincev_d_radixmerge_sort/common/include/common.hpp"
#include "votincev_d_radixmerge_sort/seq/include/ops_seq.hpp"

namespace votincev_d_radixmerge_sort {

class VotincevDRadixMergeSortRunFuncTestsProcesses : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  InType input_data;
  OutType expected_res;
  void SetUp() override {
    TestType param = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    std::string input_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_votincev_d_radixmerge_sort, param + ".txt");

    std::ifstream file(input_path);
    if (!file.is_open()) {
      return;
    }

    size_t vect_sz = 0;
    file >> vect_sz;

    std::vector<int32_t> vect_data(vect_sz);

    for (int32_t &v : vect_data) {
      file >> v;
    }

    // выставляю входные данные
    input_data = vect_data;

    // выставляю предполагаемый результат
    expected_res = vect_data;
    std::ranges::sort(expected_res);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // 1,2... процессы не владеют нужным результатом
    if (output_data.size() != expected_res.size()) {
      return true;
    }
    // 0й процесс должен иметь отсортированный массив
    return output_data == expected_res;
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

namespace {

TEST_P(VotincevDRadixMergeSortRunFuncTestsProcesses, RadixMergeSortTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {"test1", "test2", "test3", "test4", "test5",
                                             "test6", "test7", "test8", "test9", "test10"};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<VotincevDRadixMergeSortSEQ, InType>(kTestParam, PPC_SETTINGS_votincev_d_radixmerge_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    VotincevDRadixMergeSortRunFuncTestsProcesses::PrintFuncTestName<VotincevDRadixMergeSortRunFuncTestsProcesses>;

INSTANTIATE_TEST_SUITE_P(RadixSortTests, VotincevDRadixMergeSortRunFuncTestsProcesses, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace votincev_d_radixmerge_sort
