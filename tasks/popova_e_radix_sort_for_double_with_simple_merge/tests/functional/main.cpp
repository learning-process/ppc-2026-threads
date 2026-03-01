#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// #include "popova_e_radix_sort_for_double_with_simple_merge/all/include/ops_all.hpp"
#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"
// #include "popova_e_radix_sort_for_double_with_simple_merge/omp/include/ops_omp.hpp"
#include "popova_e_radix_sort_for_double_with_simple_merge/seq/include/ops_seq.hpp"
// #include "popova_e_radix_sort_for_double_with_simple_merge/stl/include/ops_stl.hpp"
// #include "popova_e_radix_sort_for_double_with_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

class PopovaERunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    int width = -1;
    int height = -1;
    int channels = -1;
    std::vector<uint8_t> img;
    // Read image in RGB to ensure consistent channel count
    {
      
    }

    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = width - height + std::min(std::accumulate(img.begin(), img.end(), 0), channels);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (input_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(PopovaERunFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

// const auto kTestTasksList =
    // std::tuple_cat(ppc::util::AddFuncTask<PopovaERadixSorForDoubleWithSimpleMergeALL, InType>(kTestParam, PPC_SETTINGS_popova_e_radix_sort_for_double_with_simple_merge),
    //                ppc::util::AddFuncTask<PopovaERadixSorForDoubleWithSimpleMergeOMP, InType>(kTestParam, PPC_SETTINGS_popova_e_radix_sort_for_double_with_simple_merge),
    //                ppc::util::AddFuncTask<PopovaERadixSorForDoubleWithSimpleMergeSEQ, InType>(kTestParam, PPC_SETTINGS_popova_e_radix_sort_for_double_with_simple_merge),
    //                ppc::util::AddFuncTask<PopovaERadixSorForDoubleWithSimpleMergeSTL, InType>(kTestParam, PPC_SETTINGS_popova_e_radix_sort_for_double_with_simple_merge),
    //                ppc::util::AddFuncTask<PopovaERadixSorForDoubleWithSimpleMergeTBB, InType>(kTestParam, PPC_SETTINGS_popova_e_radix_sort_for_double_with_simple_merge));


// const auto kTestTasksList =
//     ppc::util::AddFuncTask<PopovaERadixSorForDoubleWithSimpleMergeSEQ, InType>(kTestCases, PPC_SETTINGS_popova_e_radix_sort_for_double_with_simple_merge);
  
// const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

// const auto kPerfTestName = PopovaERunFuncTestsThreads::PrintFuncTestName<PopovaERunFuncTestsThreads>;

// INSTANTIATE_TEST_SUITE_P(RadixSortTests, PopovaERunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads