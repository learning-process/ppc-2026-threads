// tasks/peryashkin_v_binary_component_contour_processing/tests/functional/main.cpp
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"
#include "peryashkin_v_binary_component_contour_processing/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace peryashkin_v_binary_component_contour_processing {

using TestType = std::tuple<int, std::string>;

namespace {

BinaryImage MakeEmpty(int w, int h) {
  BinaryImage img;
  img.width = w;
  img.height = h;
  img.data.assign(static_cast<std::size_t>(w) * static_cast<std::size_t>(h), 0);
  return img;
}

void Set(BinaryImage &img, int x, int y, uint8_t v = 1) {
  img.data[(static_cast<std::size_t>(y) * static_cast<std::size_t>(img.width)) + static_cast<std::size_t>(x)] = v;
}

BinaryImage BuildCaseEmpty() {
  return MakeEmpty(5, 4);
}

BinaryImage BuildCasePoint() {
  auto im = MakeEmpty(5, 5);
  Set(im, 2, 2, 1);
  return im;
}

BinaryImage BuildCaseLine() {
  auto im = MakeEmpty(7, 5);
  for (int x_pos = 1; x_pos <= 5; ++x_pos) {
    Set(im, x_pos, 2, 1);
  }
  return im;
}

BinaryImage BuildCaseSquare() {
  auto im = MakeEmpty(6, 6);
  for (int y_pos = 2; y_pos <= 4; ++y_pos) {
    for (int x_pos = 2; x_pos <= 4; ++x_pos) {
      Set(im, x_pos, y_pos, 1);
    }
  }
  return im;
}

BinaryImage BuildCaseTwoComponents() {
  auto im = MakeEmpty(8, 6);
  Set(im, 1, 1);
  Set(im, 2, 1);
  Set(im, 1, 2);
  Set(im, 2, 2);
  for (int y_pos = 3; y_pos <= 5; ++y_pos) {
    Set(im, 6, y_pos);
    Set(im, 7, y_pos);
  }
  return im;
}

BinaryImage BuildCaseHole() {
  auto im = MakeEmpty(7, 7);
  for (int x_pos = 1; x_pos <= 5; ++x_pos) {
    Set(im, x_pos, 1);
    Set(im, x_pos, 5);
  }
  for (int y_pos = 1; y_pos <= 5; ++y_pos) {
    Set(im, 1, y_pos);
    Set(im, 5, y_pos);
  }
  return im;
}

BinaryImage BuildCaseTouchBorder() {
  auto im = MakeEmpty(5, 5);
  for (int y_pos = 0; y_pos <= 2; ++y_pos) {
    for (int x_pos = 0; x_pos <= 1; ++x_pos) {
      Set(im, x_pos, y_pos, 1);
    }
  }
  return im;
}

BinaryImage BuildCase(int id) {
  switch (id) {
    case 0: {
      return BuildCaseEmpty();
    }
    case 1: {
      return BuildCasePoint();
    }
    case 2: {
      return BuildCaseLine();
    }
    case 3: {
      return BuildCaseSquare();
    }
    case 4: {
      return BuildCaseTwoComponents();
    }
    case 5: {
      return BuildCaseHole();
    }
    case 6: {
      return BuildCaseTouchBorder();
    }
    default:
      return MakeEmpty(1, 1);
  }
}

}  // namespace

class PeryashkinVRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = BuildCase(std::get<0>(params));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (input_data_.data.empty()) {
      return false;
    }
    if (std::all_of(input_data_.data.begin(), input_data_.data.end(), [](uint8_t v) { return v == 0; })) {
      return output_data.empty();
    }
    return !output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
};

TEST_P(PeryashkinVRunFuncTestsThreads, BinaryComponentContourSEQ) {
  ExecuteTest(GetParam());
}

namespace {

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(0, "empty"),        std::make_tuple(1, "point"),          std::make_tuple(2, "line"),
    std::make_tuple(3, "square"),       std::make_tuple(4, "two_components"), std::make_tuple(5, "hole"),
    std::make_tuple(6, "touch_border"),
};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<PeryashkinVBinaryComponentContourProcessingSEQ, InType>(
        kTestParam, PPC_SETTINGS_peryashkin_v_binary_component_contour_processing));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kNameFn = PeryashkinVRunFuncTestsThreads::PrintFuncTestName<PeryashkinVRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(FuncTests, PeryashkinVRunFuncTestsThreads, kGtestValues, kNameFn);

}  // namespace

TEST(PeryashkinVBinaryComponentContourProcessingSEQUnit, ValidationFailsOnBadSizes) {
  BinaryImage bad;
  bad.width = 4;
  bad.height = 4;
  bad.data.assign(3, 1);
  PeryashkinVBinaryComponentContourProcessingSEQ task(bad);
  EXPECT_FALSE(task.Validation());
}

TEST(PeryashkinVBinaryComponentContourProcessingSEQUnit, ValidationFailsOnNonPositiveDims) {
  BinaryImage bad;
  bad.width = 0;
  bad.height = 5;
  bad.data.clear();
  PeryashkinVBinaryComponentContourProcessingSEQ task(bad);
  EXPECT_FALSE(task.Validation());
}

TEST(PeryashkinVBinaryComponentContourProcessingSEQUnit, PipelineReturnsFalseIfInvalid) {
  BinaryImage bad;
  bad.width = 3;
  bad.height = 3;
  bad.data.assign(1, 1);
  PeryashkinVBinaryComponentContourProcessingSEQ task(bad);
  EXPECT_FALSE(task.Validation());
  EXPECT_TRUE(task.PreProcessing());
  EXPECT_FALSE(task.Run());
}

}  // namespace peryashkin_v_binary_component_contour_processing
