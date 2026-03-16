#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "sokolov_k_matrix_double_fox/common/include/common.hpp"
#include "sokolov_k_matrix_double_fox/omp/include/ops_omp.hpp"
#include "sokolov_k_matrix_double_fox/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sokolov_k_matrix_double_fox {

class SokolovKMatrixDoubleFoxFuncTestsSeq : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(SokolovKMatrixDoubleFoxFuncTestsSeq, FoxMatmul) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {std::make_tuple(1, "single_element"),
                                             std::make_tuple(2, "2x2"),
                                             std::make_tuple(3, "3x3"),
                                             std::make_tuple(4, "4x4"),
                                             std::make_tuple(5, "5x5"),
                                             std::make_tuple(6, "6x6"),
                                             std::make_tuple(7, "7x7"),
                                             std::make_tuple(8, "8x8"),
                                             std::make_tuple(9, "9x9"),
                                             std::make_tuple(10, "10x10"),
                                             std::make_tuple(50, "50x50"),
                                             std::make_tuple(100, "100x100")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<SokolovKMatrixDoubleFoxOMP, InType>(kTestParam, PPC_SETTINGS_sokolov_k_matrix_double_fox),
    ppc::util::AddFuncTask<SokolovKMatrixDoubleFoxSEQ, InType>(kTestParam, PPC_SETTINGS_sokolov_k_matrix_double_fox));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SokolovKMatrixDoubleFoxFuncTestsSeq::PrintFuncTestName<SokolovKMatrixDoubleFoxFuncTestsSeq>;

INSTANTIATE_TEST_SUITE_P(FoxMatmulTests, SokolovKMatrixDoubleFoxFuncTestsSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sokolov_k_matrix_double_fox
