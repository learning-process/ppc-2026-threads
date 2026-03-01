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

#include "example_threads/all/include/ops_all.hpp"
#include "example_threads/common/include/common.hpp"
#include "example_threads/omp/include/ops_omp.hpp"
#include "example_threads/seq/include/ops_seq.hpp"
#include "example_threads/stl/include/ops_stl.hpp"
#include "example_threads/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace luzan_e_double_sparse_matrix_mult_seq {

class LuzanEDoubleSparseMatrixMultSeqestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    input_data_ = width - height + std::min(std::accumulate(img.begin(), img.end(), 0), channels);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(LuzanEDoubleSparseMatrixMultSeqestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<NesterovATestTaskALL, InType>(kTestParam, PPC_SETTINGS_example_threads),
    ppc::util::AddFuncTask<NesterovATestTaskOMP, InType>(kTestParam, PPC_SETTINGS_example_threads),
    ppc::util::AddFuncTask<LuzanEDoubleSparseMatrixMultSeq, InType>(kTestParam, PPC_SETTINGS_example_threads),
    ppc::util::AddFuncTask<NesterovATestTaskSTL, InType>(kTestParam, PPC_SETTINGS_example_threads),
    ppc::util::AddFuncTask<NesterovATestTaskTBB, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    LuzanEDoubleSparseMatrixMultSeqestsThreads::PrintFuncTestName<LuzanEDoubleSparseMatrixMultSeqestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, LuzanEDoubleSparseMatrixMultSeqestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace luzan_e_double_sparse_matrix_mult_seq
