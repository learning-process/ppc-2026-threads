#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
// #include <cstdint>
// #include <numeric>
// #include <stdexcept>
#include <string>
#include <tuple>
// #include <utility>
// #include <vector>

// #include "zenin_a_radix_sort_double_batcher_merge_seq/all/include/ops_all.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_seq/common/include/common.hpp"
// #include "zenin_a_radix_sort_double_batcher_merge_seq/omp/include/ops_omp.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_seq/seq/include/ops_seq.hpp"
// #include "zenin_a_radix_sort_double_batcher_merge_seq/stl/include/ops_stl.hpp"
// #include "zenin_a_radix_sort_double_batcher_merge_seq/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_seq {

class ZeninARadixSortDoubleBatcherMergeFuncTestsThreads
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType param = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(param);
    expected_data_ = std::get<1>(param);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_data_.size()) {
      return false;
    }
    for (std::size_t i = 0; i < output_data.size(); ++i) {
      if (output_data[i] != expected_data_[i]) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

namespace {

TEST_P(ZeninARadixSortDoubleBatcherMergeFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {TestType{InType{8.8}, OutType{8.8}, "SingleVal"},
                                            TestType{InType{}, OutType{}, "empty"},
                                            TestType{InType{6.6, 3.3}, OutType{3.3, 6.6}, "ReverseTwo"}};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ZeninARadixSortDoubleBatcherMergeSeqseq, InType>(
    kTestParam, PPC_SETTINGS_zenin_a_radix_sort_double_batcher_merge_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ZeninARadixSortDoubleBatcherMergeFuncTestsThreads::PrintFuncTestName<
    ZeninARadixSortDoubleBatcherMergeFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, ZeninARadixSortDoubleBatcherMergeFuncTestsThreads, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace zenin_a_radix_sort_double_batcher_merge_seq
