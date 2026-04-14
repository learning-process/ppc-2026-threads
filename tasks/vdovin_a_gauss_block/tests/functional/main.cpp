#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "vdovin_a_gauss_block/common/include/common.hpp"
#include "vdovin_a_gauss_block/omp/include/ops_omp.hpp"
#include "vdovin_a_gauss_block/seq/include/ops_seq.hpp"
#include "vdovin_a_gauss_block/tbb/include/ops_tbb.hpp"

namespace vdovin_a_gauss_block {

class VdovinAGaussBlockFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
    return output_data == 100;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(VdovinAGaussBlockFuncTestsThreads, GaussBlock) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {std::make_tuple(3, "side_3"),   std::make_tuple(4, "side_4"),
                                             std::make_tuple(5, "side_5"),   std::make_tuple(6, "side_6"),
                                             std::make_tuple(7, "side_7"),   std::make_tuple(8, "side_8"),
                                             std::make_tuple(10, "side_10"), std::make_tuple(20, "side_20"),
                                             std::make_tuple(50, "side_50"), std::make_tuple(100, "side_100")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<VdovinAGaussBlockSEQ, InType>(kTestParam, PPC_SETTINGS_vdovin_a_gauss_block),
                   ppc::util::AddFuncTask<VdovinAGaussBlockOMP, InType>(kTestParam, PPC_SETTINGS_vdovin_a_gauss_block),
                   ppc::util::AddFuncTask<VdovinAGaussBlockTBB, InType>(kTestParam, PPC_SETTINGS_vdovin_a_gauss_block));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = VdovinAGaussBlockFuncTestsThreads::PrintFuncTestName<VdovinAGaussBlockFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(GaussBlockTests, VdovinAGaussBlockFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace vdovin_a_gauss_block
