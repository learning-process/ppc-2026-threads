#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "rysev_m_linear_filter_gauss_kernel/common/include/common.hpp"
#include "rysev_m_linear_filter_gauss_kernel/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace rysev_m_linear_filter_gauss_kernel {

namespace {
OutType ComputeReferenceOutput(int test_id) {
  RysevMGaussFilterSEQ etalon(test_id);
  bool success = etalon.Validation() && etalon.PreProcessing() && etalon.Run() && etalon.PostProcessing();
  EXPECT_TRUE(success);
  return etalon.GetOutput();
}
}  // namespace

class RysevMFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);

    reference_output_ = ComputeReferenceOutput(test_id);
    input_data_ = test_id;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == reference_output_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
  OutType reference_output_ = 0;
};

namespace {

const std::array<TestType, 5> kTestParam = {std::make_tuple(0, "pic"), std::make_tuple(16, "size16"),
                                            std::make_tuple(32, "size32"), std::make_tuple(64, "size64"),
                                            std::make_tuple(128, "size128")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<RysevMGaussFilterSEQ, InType>(kTestParam, PPC_SETTINGS_rysev_m_linear_filter_gauss_kernel));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kFuncTestName = RysevMFuncTests::PrintFuncTestName<RysevMFuncTests>;

INSTANTIATE_TEST_SUITE_P(ImageTests, RysevMFuncTests, kGtestValues, kFuncTestName);

TEST_P(RysevMFuncTests, CompareWithSeq) {
  ExecuteTest(GetParam());
}

}  // namespace

}  // namespace rysev_m_linear_filter_gauss_kernel
