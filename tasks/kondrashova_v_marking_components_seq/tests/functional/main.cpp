#include <gtest/gtest.h>

#include <array>
#include <string>
#include <tuple>
#include <vector>

#include "kondrashova_v_marking_components_seq/common/include/common.hpp"
#include "kondrashova_v_marking_components_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace kondrashova_v_marking_components_seq {

class MarkingComponentsFuncTest
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType& param) {
    return std::get<1>(param);
  }

 protected:
  bool CheckTestOutputData(OutType& output_data) final {
    const std::string& type = std::get<1>(GetParam());
    
    int expected_count = 0;
    InType image = GetTestInputData();
    
    if (type == "empty") {
      expected_count = 0;
    } else if (type == "one_component") {
      expected_count = 1;
    } else if (type == "isolated_pixels") {
      expected_count = 4;
    } else if (type == "two_regions") {
      expected_count = 2;
    }

    if (output_data.count != expected_count) return false;

    if (static_cast<int>(output_data.labels.size()) != image.height) return false;
    if (!output_data.labels.empty() &&
        static_cast<int>(output_data.labels[0].size()) != image.width)
      return false;

    for (int i = 0; i < image.height; ++i) {
      for (int j = 0; j < image.width; ++j) {
        int idx = i * image.width + j;
        if (image.data[idx] == 1) {
          if (output_data.labels[i][j] != 0) return false;
        } else {
          if (output_data.labels[i][j] <= 0) return false;
        }
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    const std::string& type = std::get<1>(GetParam());
    InType image{};

    if (type == "empty") {
      image.data   = {1,1,1, 1,1,1, 1,1,1};
      image.width  = 3;
      image.height = 3;
    }
    else if (type == "one_component") {
      image.data   = {0,0,0, 0,0,0, 0,0,0};
      image.width  = 3;
      image.height = 3;
    }
    else if (type == "isolated_pixels") {
      image.data   = {0,1,0, 1,1,1, 0,1,0};
      image.width  = 3;
      image.height = 3;
    }
    else if (type == "two_regions") {
      image.data   = {
          0,0,1,1,
          0,0,1,0,
          1,1,1,0,
          1,0,0,0
      };
      image.width  = 4;
      image.height = 4;
    }

    return image;
  }
};

namespace {
TEST_P(MarkingComponentsFuncTest, VariousBinaryImages) { ExecuteTest(GetParam()); }

const std::array<TestType, 4> kTestParam = {
    std::make_tuple(0, "empty"),
    std::make_tuple(1, "one_component"),
    std::make_tuple(2, "isolated_pixels"),
    std::make_tuple(3, "two_regions")
};

const auto kTestTasksList =
    ppc::util::AddFuncTask<KondrashovaVTaskSEQ, InType>(kTestParam, PPC_SETTINGS_kondrashova_v_marking_components_seq);

INSTANTIATE_TEST_SUITE_P(MarkingComponentsFunctionalTests,
                         MarkingComponentsFuncTest,
                         ppc::util::ExpandToValues(kTestTasksList),
                         MarkingComponentsFuncTest::PrintFuncTestName<MarkingComponentsFuncTest>);
}  // namespace
}  // namespace kondrashova_v_marking_components_seq