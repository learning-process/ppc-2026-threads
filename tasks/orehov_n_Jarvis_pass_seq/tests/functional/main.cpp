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

#include "orehov_n_Jarvis_pass_seq/common/include/common.hpp"
#include "orehov_n_Jarvis_pass_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace orehov_n_Jarvis_pass_seq {

class OrehovNJarvisPassSEQFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param);
  }

 protected:
  void SetUp() override {
    TestType param = std::get<static_cast<int>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    if (param == 1){
      input_data_.push_back(Point(0, 0));
      input_data_.push_back(Point(1, 0));
      input_data_.push_back(Point(3, 0));
      input_data_.push_back(Point(2, 1));
      input_data_.push_back(Point(4, 2));
      input_data_.push_back(Point(0, 3));
      input_data_.push_back(Point(2, 4));
      _test_res.push_back(Point(0, 0));
      _test_res.push_back(Point(0, 3));
      _test_res.push_back(Point(2, 4));
      _test_res.push_back(Point(4, 2));
      _test_res.push_back(Point(3, 0));
    }
    if (param == 2){
      input_data_.push_back(Point(0, 0));
      _test_res.push_back(Point(0, 0));
    }
    if (param == 3){
      input_data_.push_back(Point(0, 0));
      input_data_.push_back(Point(1, 0));
      _test_res.push_back(Point(0, 0));
      _test_res.push_back(Point(1, 0));
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (_test_res == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = std::vector<Point>();
  std::vector<Point> _test_res;
};

namespace {

TEST_P(OrehovNJarvisPassSEQFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {1, 2, 3};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<OrehovNJarvisPassSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = OrehovNJarvisPassSEQFuncTests::PrintFuncTestName<OrehovNJarvisPassSEQFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, OrehovNJarvisPassSEQFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace orehov_n_Jarvis_pass_seq
