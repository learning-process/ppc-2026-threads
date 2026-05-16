#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <fstream>
#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include "lopatin_a_sobel_operator/common/include/common.hpp"
#include "lopatin_a_sobel_operator/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace lopatin_a_sobel_operator {

class LopatinARunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string filename = params + ".txt";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_lopatin_a_sobel_operator, filename);
    std::ifstream infile(abs_path, std::ios::in);
    if (!infile.is_open()) {
      throw std::runtime_error("Failed to open file: " + filename);
    }

    std::string line;
    std::size_t count = 0;
    while (std::getline(infile, line)) {
      if (line.empty()) {
        continue;
      }

      std::istringstream iss(line);
      int elem = 0;

      if (count == 0) {
        iss >> input_data_.height;
      } else if (count == 1) {
        iss >> input_data_.width;
      } else if (count < input_data_.height + 2) {
        while (iss >> elem) {
          input_data_.pixels.push_back(elem);
        }
      } else {
        while (iss >> elem) {
          output_chekup_data_.push_back(elem);
        }
      }

      ++count;
    }

    input_data_.threshold = 100;

    infile.close();
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == output_chekup_data_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType output_chekup_data_;
};

namespace {

TEST_P(LopatinARunFuncTests, SobelOperator) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::string("test_img_3x3"), std::string("test_img_5x5"),
                                            std::string("test_img_7x7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<LopatinASobelOperatorSEQ, InType>(kTestParam, PPC_SETTINGS_lopatin_a_sobel_operator));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = LopatinARunFuncTests::PrintFuncTestName<LopatinARunFuncTests>;

INSTANTIATE_TEST_SUITE_P(SobelOperatorTests, LopatinARunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace lopatin_a_sobel_operator
