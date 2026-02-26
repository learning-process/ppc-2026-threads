#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"
#include "frolova_s_radix_sort_double/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace frolova_s_radix_sort_double {

class FrolovaSRadixSortDoubleRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  InType input_data;
  OutType expected_res;

  void SetUp() override {
    TestType param = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    std::string input_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_frolova_s_radix_sort_double, param + ".txt");

    std::ifstream file(input_path);
    if (!file.is_open()) {
      return;
    }

    size_t vect_sz = 0;
    file >> vect_sz;

    std::vector<double> vect_data(vect_sz);
    for (double &v : vect_data) {
      file >> v;
    }

    input_data = vect_data;
    expected_res = vect_data;
    std::sort(expected_res.begin(), expected_res.end());
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Для последовательной версии размер должен совпадать
    if (output_data.size() != expected_res.size()) {
      return false;
    }
    return output_data == expected_res;
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

namespace {

TEST_P(FrolovaSRadixSortDoubleRunFuncTests, RadixSortDoubleTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {
    "test1", "test2", "test3", "test4", "test5",
    "test6", "test7", "test8", "test9", "test10"
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<FrolovaSRadixSortDoubleSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = FrolovaSRadixSortDoubleRunFuncTests::PrintFuncTestName<FrolovaSRadixSortDoubleRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(RadixSortDoubleTests, FrolovaSRadixSortDoubleRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace frolova_s_radix_sort_double