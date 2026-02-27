#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// #include "golovanov_d_radix_merge/all/include/ops_all.hpp"
#include "golovanov_d_radix_merge/common/include/common.hpp"
// #include "golovanov_d_radix_merge/omp/include/ops_omp.hpp"
#include "golovanov_d_radix_merge/seq/include/ops_seq.hpp"
// #include "golovanov_d_radix_merge/stl/include/ops_stl.hpp"
// #include "golovanov_d_radix_merge/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace golovanov_d_radix_merge {

class GolovanovDRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
  InType input_data_{};

 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param.size());
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = params;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

namespace {

TEST_P(GolovanovDRunFuncTestsThreads, RadixMergeFunc) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {TestType{
                                                3.14,
                                                -2.71,
                                                0.0,
                                                42.0,
                                                -1.5,
                                            },
                                            TestType{-10.0, -5.0, -3.0, -1.0, -7.0, 6},
                                            TestType{100.5, 1.1, 50.0, 2.2, 0.0001, 12, 2}};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<GolovanovDRadixMergeSEQ, InType>(kTestParam, PPC_SETTINGS_golovanov_d_radix_merge));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = GolovanovDRunFuncTestsThreads::PrintFuncTestName<GolovanovDRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(RadixMergeFunc, GolovanovDRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace golovanov_d_radix_merge
