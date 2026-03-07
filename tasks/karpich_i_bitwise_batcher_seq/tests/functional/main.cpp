#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "karpich_i_bitwise_batcher_seq/common/include/common.hpp"
#include "karpich_i_bitwise_batcher_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace karpich_i_bitwise_batcher_seq {

class KarpichIBitwiseBatcherFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
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

TEST_P(KarpichIBitwiseBatcherFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {
    std::make_tuple(3, "3"),     std::make_tuple(5, "5"),       std::make_tuple(7, "7"),
    std::make_tuple(1, "1"),     std::make_tuple(2, "2"),       std::make_tuple(4, "4"),
    std::make_tuple(8, "8"),     std::make_tuple(16, "16"),     std::make_tuple(100, "100"),
    std::make_tuple(256, "256"), std::make_tuple(1000, "1000"), std::make_tuple(1024, "1024"),
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KarpichIBitwiseBatcherSEQ, InType>(kTestParam, PPC_SETTINGS_karpich_i_bitwise_batcher_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    KarpichIBitwiseBatcherFuncTestsThreads::PrintFuncTestName<KarpichIBitwiseBatcherFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, KarpichIBitwiseBatcherFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace karpich_i_bitwise_batcher_seq
