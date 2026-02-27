#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"
#include "papulina_y_radix_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace papulina_y_radix_sort {

class PapulinaYRunFuncTestThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_result_ = std::get<0>(params);
    std::ranges::sort(expected_result_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data == expected_result_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  std::vector<double> expected_result_;
};

namespace {

TEST_P(PapulinaYRunFuncTestThreads, RadixSortTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(std::vector<double>{5.0, 1.2, 0.0, 3.0, 2.0, 0.2}, "test1"),
    std::make_tuple(std::vector<double>{-5.0, -1.2, 0.0, -3.0, -2.0, -0.2}, "test2"),
    std::make_tuple(std::vector<double>{5.0, -1.2, 0.0, 3.0, -2.0, 0.2}, "test3"),
    std::make_tuple(std::vector<double>{10.0, 9.1234, -3.45667, -5.6578, 8.13154, 7.09876, 0.0, 0.12454, -0.1232},
                    "test4"),
    std::make_tuple(std::vector<double>{5.55555, 5.55554, 5.55553, 5.555545, 5.555551, 5.555556}, "test5"),
    std::make_tuple(
        std::vector<double>{
            0.00000001,
            0.000000011,
            0.000000011,
            0.0000000112,
            0.00000001112,
            0.00000001111,
            0.00000001111,
        },
        "test6"),
    std::make_tuple(std::vector<double>(10, 0.123456789), "test7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<PapulinaYRadixSortSEQ, InType>(kTestParam, PPC_SETTINGS_papulina_y_radix_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = PapulinaYRunFuncTestThreads::PrintFuncTestName<PapulinaYRunFuncTestThreads>;

INSTANTIATE_TEST_SUITE_P(RadixSortTests, PapulinaYRunFuncTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace papulina_y_radix_sort
