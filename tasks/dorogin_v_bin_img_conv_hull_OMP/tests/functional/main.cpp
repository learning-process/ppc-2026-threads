#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_OMP/common/include/common.hpp"
#include "dorogin_v_bin_img_conv_hull_OMP/omp/include/ops_omp.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dorogin_v_bin_img_conv_hull_omp {

namespace {

BinaryImage MakeEmpty(int w, int h) {
  BinaryImage img;
  img.width = w;
  img.height = h;
  img.data.assign(static_cast<std::size_t>(w) * static_cast<std::size_t>(h), 0);
  return img;
}

void SetPixel(BinaryImage &img, int x, int y, std::uint8_t v = 1U) {
  img.data[(static_cast<std::size_t>(y) * static_cast<std::size_t>(img.width)) + static_cast<std::size_t>(x)] = v;
}

BinaryImage CaseSinglePoint() {
  auto img = MakeEmpty(5, 5);
  SetPixel(img, 2, 2, 1U);
  return img;
}

BinaryImage CaseHorizontalSegment() {
  auto img = MakeEmpty(7, 3);
  for (int xx = 1; xx <= 5; ++xx) {
    SetPixel(img, xx, 1, 1U);
  }
  return img;
}

BinaryImage CaseRectangle() {
  auto img = MakeEmpty(6, 6);
  for (int yy = 1; yy <= 4; ++yy) {
    for (int xx = 2; xx <= 4; ++xx) {
      SetPixel(img, xx, yy, 1U);
    }
  }
  return img;
}

BinaryImage CaseTwoSeparated() {
  auto img = MakeEmpty(8, 5);
  SetPixel(img, 1, 1);
  SetPixel(img, 2, 1);
  SetPixel(img, 1, 2);
  SetPixel(img, 2, 2);
  for (int yy = 0; yy < 5; ++yy) {
    SetPixel(img, 6, yy);
  }
  return img;
}

BinaryImage CaseAllBackground() {
  return MakeEmpty(4, 4);
}

BinaryImage BuildCase(int id) {
  switch (id) {
    case 0:
      return CaseAllBackground();
    case 1:
      return CaseSinglePoint();
    case 2:
      return CaseHorizontalSegment();
    case 3:
      return CaseRectangle();
    case 4:
      return CaseTwoSeparated();
    default:
      return MakeEmpty(1, 1);
  }
}

}  // namespace

class DoroginVRunFuncTestsOMP : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
    const bool has_foreground = std::ranges::any_of(input_data_.data, [](std::uint8_t v) { return v != 0U; });
    if (!has_foreground) {
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

namespace {

TEST_P(DoroginVRunFuncTestsOMP, BinaryImageConvexHullsOMP) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(0, "background"), std::make_tuple(1, "single_point"),   std::make_tuple(2, "segment"),
    std::make_tuple(3, "rectangle"),  std::make_tuple(4, "two_components"),
};

const auto kOmpTasksList = std::tuple_cat(ppc::util::AddFuncTask<DoroginVBinImgConvHullOMP, InType>(
    kTestParam, PPC_SETTINGS_dorogin_v_bin_img_conv_hull_OMP));

const auto kGtestValues = ppc::util::ExpandToValues(kOmpTasksList);

const auto kTestName = DoroginVRunFuncTestsOMP::PrintFuncTestName<DoroginVRunFuncTestsOMP>;

INSTANTIATE_TEST_SUITE_P(FuncTests, DoroginVRunFuncTestsOMP, kGtestValues, kTestName);

}  // namespace

}  // namespace dorogin_v_bin_img_conv_hull_omp
