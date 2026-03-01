#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"
#include "moskaev_v_lin_filt_block_gauss_3/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

class MoskaevVLinFiltBlockGauss3SEQFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<6>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::make_tuple(std::get<0>(params), std::get<1>(params), std::get<2>(params), std::get<3>(params),
                                  std::get<4>(params));
    expected_result_ = std::get<5>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_result_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_result_;
};

namespace {

TEST_P(MoskaevVLinFiltBlockGauss3SEQFuncTests, GaussFilter) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {
    std::make_tuple(2, 2, 1, 64, std::vector<uint8_t>{100, 150, 200, 250}, std::vector<uint8_t>{138, 163, 188, 213},
                    std::string("test1_2x2_gray")),

    std::make_tuple(3, 3, 1, 64, std::vector<uint8_t>{1, 2, 3, 4, 5, 6, 7, 8, 9},
                    std::vector<uint8_t>{2, 3, 4, 4, 5, 6, 7, 7, 8}, std::string("test2_3x3_gray")),

    std::make_tuple(2, 2, 3, 64, std::vector<uint8_t>{10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120},
                    std::vector<uint8_t>{33, 43, 53, 48, 58, 68, 63, 73, 83, 78, 88, 98}, std::string("test3_2x2_rgb")),

    std::make_tuple(1, 1, 1, 64, std::vector<uint8_t>{255}, std::vector<uint8_t>{255}, std::string("test4_1x1"))};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<MoskaevVLinFiltBlockGauss3SEQ, InType>(
    kTestParam, PPC_SETTINGS_moskaev_v_lin_filt_block_gauss_3));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName =
    MoskaevVLinFiltBlockGauss3SEQFuncTests::PrintFuncTestName<MoskaevVLinFiltBlockGauss3SEQFuncTests>;

INSTANTIATE_TEST_SUITE_P(GaussFilter, MoskaevVLinFiltBlockGauss3SEQFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace moskaev_v_lin_filt_block_gauss_3
