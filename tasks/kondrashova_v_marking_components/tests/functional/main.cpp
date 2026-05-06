#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "kondrashova_v_marking_components/common/include/common.hpp"
#include "kondrashova_v_marking_components/omp/include/ops_omp.hpp"
#include "kondrashova_v_marking_components/seq/include/ops_seq.hpp"
#include "kondrashova_v_marking_components/stl/include/ops_stl.hpp"
#include "kondrashova_v_marking_components/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"

namespace kondrashova_v_marking_components {

namespace {

int GetExpectedCount(const std::string &type) {
  if (type == "zero_width" || type == "zero_height") {
    return 0;
  }
  if (type == "one_component") {
    return 1;
  }
  if (type == "isolated_pixels") {
    return 4;
  }
  if (type == "two_regions") {
    return 2;
  }
  if (type == "u_shape") {
    return 1;
  }
  if (type == "complex") {
    return 2;
  }
  if (type == "single_row_gaps" || type == "single_column_gaps") {
    return 2;
  }
  if (type == "merge_labels" || type == "boundary_bridge") {
    return 1;
  }
  if (type == "large_complex") {
    return 4;
  }
  return 0;
}

bool CheckLabelsSize(const OutType &output_data, const InType &image) {
  if (image.width <= 0 || image.height <= 0) {
    return output_data.labels.empty();
  }
  if (output_data.labels.size() != static_cast<size_t>(image.height)) {
    return false;
  }
  if (!output_data.labels.empty() && output_data.labels[0].size() != static_cast<size_t>(image.width)) {
    return false;
  }
  return true;
}

bool CheckLabelsValues(const OutType &output_data, const InType &image) {
  if (image.width <= 0 || image.height <= 0) {
    return true;
  }
  for (int ii = 0; ii < image.height; ++ii) {
    for (int jj = 0; jj < image.width; ++jj) {
      auto idx = (static_cast<size_t>(ii) * static_cast<size_t>(image.width)) + static_cast<size_t>(jj);
      const bool is_background = image.data[idx] == 1;
      const bool label_is_valid = is_background ? (output_data.labels[ii][jj] == 0) : (output_data.labels[ii][jj] > 0);
      if (!label_is_valid) {
        return false;
      }
    }
  }
  return true;
}

// 4 components to defeat GCC -O3 unrolling (kept out of GetTestInputData for clang-tidy complexity).
InType MakeLargeComplexTestImage() {
  constexpr int kSize = 30;
  InType image{};
  image.width = kSize;
  image.height = kSize;
  image.data.assign(static_cast<size_t>(kSize) * static_cast<size_t>(kSize), 1);

  for (int i = 0; i < kSize; ++i) {
    image.data[static_cast<size_t>(i)] = 0;
  }

  for (int i = 5; i < 25; ++i) {
    image.data[(i * kSize) + 15] = 0;
  }

  for (int i = 25; i <= 27; ++i) {
    for (int j = 25; j <= 27; ++j) {
      image.data[(i * kSize) + j] = 0;
    }
  }

  for (int i = 10; i < 20; ++i) {
    image.data[(i * kSize) + 5] = 0;  // left arm
  }
  for (int i = 10; i < 20; ++i) {
    image.data[(i * kSize) + 10] = 0;  // right arm
  }
  for (int j = 5; j <= 10; ++j) {
    image.data[(19 * kSize) + j] = 0;  // bottom connection
  }

  return image;
}

}  // namespace

class MarkingComponentsFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    return std::get<1>(param);
  }

 protected:
  bool CheckTestOutputData(OutType &output_data) final {
    const std::string &type = std::get<1>(GetParam());
    const InType image = GetTestInputData();

    if (output_data.count != GetExpectedCount(type)) {
      return false;
    }
    if (!CheckLabelsSize(output_data, image)) {
      return false;
    }
    if (!CheckLabelsValues(output_data, image)) {
      return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    const std::string &type = std::get<1>(GetParam());
    InType image{};

    if (type == "zero_width") {
      image.width = 0;
      image.height = 3;
    } else if (type == "zero_height") {
      image.width = 3;
      image.height = 0;
    } else if (type == "empty") {
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
    } else if (type == "u_shape") {
      image.data = {0, 1, 0, 0, 1, 0, 0, 0, 0};
      image.width = 3;
      image.height = 3;
    } else if (type == "complex") {
      // 2 components
      image.data = {0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1};
      image.width = 5;
      image.height = 5;
    } else if (type == "single_row_gaps") {
      image.data = {0, 0, 1, 0, 0};
      image.width = 5;
      image.height = 1;
    } else if (type == "single_column_gaps") {
      image.data = {0, 0, 1, 0, 0};
      image.width = 1;
      image.height = 5;
    } else if (type == "merge_labels") {
      image.data = {0, 1, 0, 0, 0, 0, 1, 1, 1};
      image.width = 3;
      image.height = 3;
    } else if (type == "boundary_bridge") {
      image.data = {1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1};
      image.width = 3;
      image.height = 8;
    } else if (type == "large_complex") {
      image = MakeLargeComplexTestImage();
    }

    return image;
  }
};

namespace {
TEST_P(MarkingComponentsFuncTest, VariousBinaryImages) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 13> kTestParam = {
    std::make_tuple(0, "zero_width"),      std::make_tuple(1, "zero_height"),
    std::make_tuple(2, "empty"),           std::make_tuple(3, "one_component"),
    std::make_tuple(4, "isolated_pixels"), std::make_tuple(5, "two_regions"),
    std::make_tuple(6, "u_shape"),         std::make_tuple(7, "complex"),
    std::make_tuple(8, "single_row_gaps"), std::make_tuple(9, "single_column_gaps"),
    std::make_tuple(10, "merge_labels"),   std::make_tuple(11, "boundary_bridge"),
    std::make_tuple(12, "large_complex")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KondrashovaVTaskSEQ, InType>(kTestParam, PPC_SETTINGS_kondrashova_v_marking_components),
    ppc::util::AddFuncTask<KondrashovaVTaskOMP, InType>(kTestParam, PPC_SETTINGS_kondrashova_v_marking_components),
    ppc::util::AddFuncTask<KondrashovaVTaskSTL, InType>(kTestParam, PPC_SETTINGS_kondrashova_v_marking_components),
    ppc::util::AddFuncTask<KondrashovaVTaskTBB, InType>(kTestParam, PPC_SETTINGS_kondrashova_v_marking_components));

INSTANTIATE_TEST_SUITE_P(KondrashovaVMarkingComponentsFunctionalTests, MarkingComponentsFuncTest,
                         ppc::util::ExpandToValues(kTestTasksList),
                         MarkingComponentsFuncTest::PrintFuncTestName<MarkingComponentsFuncTest>);
}  // namespace

}  // namespace kondrashova_v_marking_components
