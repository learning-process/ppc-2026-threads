#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <tuple>

#include "egashin_k_radix_simple_merge/all/include/ops_all.hpp"
#include "egashin_k_radix_simple_merge/common/include/common.hpp"
#include "egashin_k_radix_simple_merge/omp/include/ops_omp.hpp"
#include "egashin_k_radix_simple_merge/seq/include/ops_seq.hpp"
#include "egashin_k_radix_simple_merge/stl/include/ops_stl.hpp"
#include "egashin_k_radix_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace egashin_k_radix_simple_merge {

class EgashinKRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_result_ = input_data_;
    std::ranges::sort(expected_result_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_result_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_result_;
};

namespace {

TEST_P(EgashinKRunFuncTestsThreads, SortsDoubleVector) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {{
    {InType{}, "empty_vector"},
    {InType{42.0}, "single_element"},
    {InType{-10.5, -1.0, 0.25, 5.0, 12.75}, "already_sorted"},
    {InType{9.0, 7.0, 5.0, 3.0, 1.0}, "reverse_sorted"},
    {InType{5.0, -1.25, 0.0, 3.5, -2.75, 0.25}, "mixed_signs"},
    {InType{4.5, 4.5, -2.0, 4.5, -2.0, 8.1, 8.1}, "duplicates"},
    {InType{1e-12, -1e12, std::numbers::pi, -std::numbers::e, 1e6, -1e-6}, "different_scales"},
    {InType{std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -10.0, 1.0, 0.5},
     "with_infinities"},
}};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<EgashinKRadixSimpleMergeSEQ, InType>(kTestParam, PPC_SETTINGS_egashin_k_radix_simple_merge),
    ppc::util::AddFuncTask<EgashinKRadixSimpleMergeOMP, InType>(kTestParam, PPC_SETTINGS_egashin_k_radix_simple_merge),
    ppc::util::AddFuncTask<EgashinKRadixSimpleMergeTBB, InType>(kTestParam, PPC_SETTINGS_egashin_k_radix_simple_merge),
    ppc::util::AddFuncTask<EgashinKRadixSimpleMergeSTL, InType>(kTestParam, PPC_SETTINGS_egashin_k_radix_simple_merge),
    ppc::util::AddFuncTask<EgashinKRadixSimpleMergeALL, InType>(kTestParam, PPC_SETTINGS_egashin_k_radix_simple_merge));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = EgashinKRunFuncTestsThreads::PrintFuncTestName<EgashinKRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(RadixSortTests, EgashinKRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace egashin_k_radix_simple_merge
