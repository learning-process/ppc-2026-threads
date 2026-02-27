#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "kondrashova_v_marking_components_seq/common/include/common.hpp"
#include "kondrashova_v_marking_components_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace kondrashova_v_marking_components_seq {

namespace {

int GetExpectedCount(const std::string &type) {
  if (type == "one_component") return 1;
  if (type == "isolated_pixels") return 4;
  if (type == "two_regions") return 2;
  return 0;
}

bool CheckLabelsSize(const OutType &output_data, const InType &image) {
  if (output_data.labels.size() != static_cast<size_t>(image.height)) {
    return false;
  }
  if (!output_data.labels.empty() &&
      output_data.labels[0].size() != static_cast<size_t>(image.width)) {
    return false;
  }
  return true;
}

bool CheckLabelsValues(const OutType &output_data, const InType &image) {
  for (int ii = 0; ii < image.height; ++ii) {
    for (int jj = 0; jj < image.width; ++jj) {
      auto idx = static_cast<size_t>(ii) * static_cast<size_t>(image.width) +
                 static_cast<size_t>(jj);
      if (image.data[idx] == 1) {
        if (output_data.labels[ii][jj] != 0) return false;
      } else {
        if (output_data.labels[ii][jj] <= 0) return false;
      }
    }
  }
  return true;
}

}  // namespace

class MarkingComponentsFuncTest
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    return std::get<1>(param);
  }

 protected:
  bool CheckTestOutputData(OutType &output_data) final {
    const std::string &type = std::get<1>(GetParam());
    const InType image = GetTestInputData();

    if (output_data.count != GetExpectedCount(type)) return false;
    if (!CheckLabelsSize(output_data, image)) return false;
    if (!CheckLabelsValues(output_data, image)) return false;
    return true;
  }

  InType GetTestInputData() final {
    const std::string &type = std::get<1>(GetParam());
    InType image{};

    if (type == "empty") {
      image.data = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      image.width = 3;
      image.height = 3;
    } else if (type == "one_component") {
      image.data = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      image.width = 3;
      image.height = 3;
    } else if (type == "isolated_pixels") {
      image.data = {0, 1, 0, 1, 1, 1, 0, 1, 0};
      image.width = 3;
      image.height = 3;
    } else if (type == "two_regions") {
      image.data = {0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0};
      image.width = 4;
      image.height = 4;
    }

    return image;
  }
};

namespace {
TEST_P(MarkingComponentsFuncTest, VariousBinaryImages) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {
    std::make_tuple(0, "empty"), std::make_tuple(1, "one_component"),
    std::make_tuple(2, "isolated_pixels"), std::make_tuple(3, "two_regions")};

const auto kTestTasksList = ppc::util::AddFuncTask<KondrashovaVTaskSEQ, InType>(
    kTestParam, PPC_SETTINGS_kondrashova_v_marking_components_seq);

INSTANTIATE_TEST_SUITE_P(
    KondrashovaVMarkingComponentsFunctionalTests, MarkingComponentsFuncTest,
    ppc::util::ExpandToValues(kTestTasksList),
    MarkingComponentsFuncTest::PrintFuncTestName<MarkingComponentsFuncTest>);
}  // namespace

}  // namespace kondrashova_v_marking_components_seq