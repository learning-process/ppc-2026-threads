#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
#include "gonozov_l_bitwise_sorting_double_Batcher_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace gonozov_l_bitwise_sorting_double_Batcher_merge {

class GonozovLBitSortBatcherMergeFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  InType input_data_;

  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::vector<double> ordered_data = static_cast<OutType>(std::get<1>(params));
    size_t num_elem = ordered_data.size();
    for (size_t i = 0; i < num_elem; i++) {
      if (std::abs(ordered_data[i] - output_data[i]) > 0.01) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

namespace {

TEST_P(GonozovLBitSortBatcherMergeFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {
    std::make_tuple(std::vector<double>{3.0, 7.0, 1.0, 5.0, 2.0, 11.0},
                    std::vector<double>{1.0, 2.0, 3.0, 5.0, 7.0, 11.0}, "sorting_positive_simple_integers"),
    std::make_tuple(std::vector<double>{3.5, 7.6, 7.1, 5.8, 5.2, 1.1},
                    std::vector<double>{1.1, 3.5, 5.2, 5.8, 7.1, 7.6}, "sorting_positive_double"),
    std::make_tuple(std::vector<double>{-3.0, 7.0, -1.0, 5.0, -2.0, 11.0},
                    std::vector<double>{-3.0, -2.0, -1.0, 5.0, 7.0, 11.0}, "sorting_pos_and_neg_integers"),
    std::make_tuple(std::vector<double>{3.5, -7.6, 7.1, -5.8, 5.2, -1.1},
                    std::vector<double>{-7.6, -5.8, -1.1, 3.5, 5.2, 7.1}, "sorting_pos_and_neg_double"),
};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<GonozovLBitSortBatcherMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_gonozov_l_bitwise_sorting_double_Batcher_merge));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    GonozovLBitSortBatcherMergeFuncTests::PrintFuncTestName<GonozovLBitSortBatcherMergeFuncTests>;

INSTANTIATE_TEST_SUITE_P(FuncTests, GonozovLBitSortBatcherMergeFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace gonozov_l_bitwise_sorting_double_Batcher_merge
