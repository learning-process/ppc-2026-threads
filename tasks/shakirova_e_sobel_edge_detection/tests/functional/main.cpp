#include <gtest/gtest.h>

#include <array>
#include <cctype>
#include <cstddef>
#include <string>
#include <tuple>

#include "shakirova_e_sobel_edge_detection/common/include/common.hpp"
#include "shakirova_e_sobel_edge_detection/common/include/img_container.hpp"
#include "shakirova_e_sobel_edge_detection/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shakirova_e_sobel_edge_detection {

class ShakirovaESobelEdgeDetectionFuncTestsThreads
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType& test_param) {
    std::string name = std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
    for (char& c : name) {
      if ((std::isalnum(static_cast<unsigned char>(c)) == 0) && c != '_') {
        c = '_';
      }
    }
    return name;
  }

 protected:
  void SetUp() override {
    const auto& [expected, filename] =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    expected_edges_ = expected;

    const std::string abs_path =
        ppc::util::GetAbsoluteTaskPath(
            std::string(PPC_ID_shakirova_e_sobel_edge_detection),
            filename);
    input_data_ = ImgContainer::FromFile(abs_path);
  }

  bool CheckTestOutputData(OutType& output_data) final {
    int edge_count = 0;
    for (int row = 1; row < input_data_.height - 1; ++row) {
      for (int col = 1; col < input_data_.width - 1; ++col) {
        if (output_data[(row * input_data_.width) + col] > 0) {
          ++edge_count;
        }
      }
    }
    return (expected_edges_ == 0) ? (edge_count == 0) : (edge_count >= expected_edges_);
  }

  InType GetTestInputData() final { return input_data_; }

 private:
  InType input_data_;
  int expected_edges_{0};
};

namespace {

TEST_P(ShakirovaESobelEdgeDetectionFuncTestsThreads, SobelOnFiles) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParam = {
    std::make_tuple(4,  "test_3.png"),
    std::make_tuple(10, "test_1.png"),
    std::make_tuple(1,  "test_2.png"),
    std::make_tuple(0,  "test_4.txt"),
    std::make_tuple(2,  "test_5.txt"),
    std::make_tuple(2,  "test_6.txt"),
};

const auto kTestTasksList =
    ppc::util::AddFuncTask<ShakirovaESobelEdgeDetectionSEQ, InType>(
        kTestParam, PPC_SETTINGS_shakirova_e_sobel_edge_detection);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    ShakirovaESobelEdgeDetectionFuncTestsThreads::PrintFuncTestName<ShakirovaESobelEdgeDetectionFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SobelEdgeTests, ShakirovaESobelEdgeDetectionFuncTestsThreads,
                         kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shakirova_e_sobel_edge_detection