#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/common/include/test_reader.hpp"
#include "guseva_crs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace guseva_crs {

class GusevaMatMulCRSFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    const auto &filename = ppc::util::GetAbsoluteTaskPath(
        std::string(PPC_ID_guseva_crs),
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()) + ".txt");
    const auto &[a, b, c] = ReadTestFromFile(filename);
    input_data_ = std::make_tuple(a, b);
    output_data_ = c;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return Equal(output_data_, output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType output_data_;
};

namespace {

TEST_P(GusevaMatMulCRSFuncTest, G) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParam = {"sparse_dense",   "dense_sparse",   "double_sparse1",
                                            "double_sparse2", "double_sparse3", "double_sparse4"};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<GusevaCRSMatMulSeq, InType>(kTestParam, PPC_SETTINGS_guseva_crs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = GusevaMatMulCRSFuncTest::PrintFuncTestName<GusevaMatMulCRSFuncTest>;

INSTANTIATE_TEST_SUITE_P(MatMulSRC, GusevaMatMulCRSFuncTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace guseva_crs
