#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm> // abc
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "rysev_m_linear_filter_gauss_kernel/all/include/ops_all.hpp"
#include "rysev_m_linear_filter_gauss_kernel/common/include/common.hpp"
#include "rysev_m_linear_filter_gauss_kernel/omp/include/ops_omp.hpp"
#include "rysev_m_linear_filter_gauss_kernel/seq/include/ops_seq.hpp"
#include "rysev_m_linear_filter_gauss_kernel/stl/include/ops_stl.hpp"
#include "rysev_m_linear_filter_gauss_kernel/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace rysev_m_linear_filter_gauss_kernel {

class RysevMFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    RysevMGaussFilterSEQ etalon(0);
    ASSERT_TRUE(etalon.Validation());
    ASSERT_TRUE(etalon.PreProcessing());
    ASSERT_TRUE(etalon.Run());
    ASSERT_TRUE(etalon.PostProcessing());
    reference_output_ = etalon.GetOutput();

    input_data_ = 0;
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

TEST_P(RysevMFuncTests, CompareWithSeq) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 1> kTestParam = {std::make_tuple(0, "pic")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<RysevMGaussFilterALL, InType>(kTestParam, PPC_SETTINGS_rysev_m_linear_filter_gauss_kernel),
                   ppc::util::AddFuncTask<RysevMGaussFilterOMP, InType>(kTestParam, PPC_SETTINGS_rysev_m_linear_filter_gauss_kernel),
                   ppc::util::AddFuncTask<RysevMGaussFilterSEQ, InType>(kTestParam, PPC_SETTINGS_rysev_m_linear_filter_gauss_kernel),
                   ppc::util::AddFuncTask<RysevMGaussFilterSTL, InType>(kTestParam, PPC_SETTINGS_rysev_m_linear_filter_gauss_kernel),
                   ppc::util::AddFuncTask<RysevMGaussFilterTBB, InType>(kTestParam, PPC_SETTINGS_rysev_m_linear_filter_gauss_kernel));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = RysevMFuncTests::PrintFuncTestName<RysevMFuncTests>;

INSTANTIATE_TEST_SUITE_P(ImageTests, RysevMFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace rysev_m_linear_filter_gauss_kernel