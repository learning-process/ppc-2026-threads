#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"
#include "peryashkin_v_binary_component_contour_processing/omp/include/ops_omp.hpp"
#include "peryashkin_v_binary_component_contour_processing/seq/include/ops_seq.hpp"
#include "peryashkin_v_binary_component_contour_processing/stl/include/ops_stl.hpp"
#include "peryashkin_v_binary_component_contour_processing/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace peryashkin_v_binary_component_contour_processing {

using TestType = std::tuple<int, std::string>;

namespace {

BinaryImage MakeEmpty(int w, int h) {
  BinaryImage img;
  img.width = w;
  img.height = h;
  img.data.assign((static_cast<std::size_t>(w) * static_cast<std::size_t>(h)), 0);
  return img;
}

void Set(BinaryImage &img, int x, int y, std::uint8_t v = 1) {
  img.data[(static_cast<std::size_t>(y) * static_cast<std::size_t>(img.width)) + static_cast<std::size_t>(x)] = v;
}

BinaryImage CaseEmpty() {
  return MakeEmpty(5, 4);
}

BinaryImage CasePoint() {
  auto im = MakeEmpty(5, 5);
  Set(im, 2, 2, 1);
  return im;
}

BinaryImage CaseLine() {
  auto im = MakeEmpty(7, 5);
  for (int xx = 1; xx <= 5; ++xx) {
    Set(im, xx, 2, 1);
  }
  return im;
}

BinaryImage CaseSquare() {
  auto im = MakeEmpty(6, 6);
  for (int yy = 2; yy <= 4; ++yy) {
    for (int xx = 2; xx <= 4; ++xx) {
      Set(im, xx, yy, 1);
    }
  }
  return im;
}

BinaryImage CaseTwoComponents() {
  auto im = MakeEmpty(8, 6);
  Set(im, 1, 1);
  Set(im, 2, 1);
  Set(im, 1, 2);
  Set(im, 2, 2);
  for (int yy = 3; yy <= 5; ++yy) {
    Set(im, 6, yy);
    Set(im, 7, yy);
  }
  return im;
}

BinaryImage CaseHole() {
  auto im = MakeEmpty(7, 7);
  for (int xx = 1; xx <= 5; ++xx) {
    Set(im, xx, 1);
    Set(im, xx, 5);
  }
  for (int yy = 1; yy <= 5; ++yy) {
    Set(im, 1, yy);
    Set(im, 5, yy);
  }
  return im;
}

BinaryImage CaseTouchBorder() {
  auto im = MakeEmpty(5, 5);
  for (int yy = 0; yy <= 2; ++yy) {
    for (int xx = 0; xx <= 1; ++xx) {
      Set(im, xx, yy, 1);
    }
  }
  return im;
}

BinaryImage BuildCase(int id) {
  switch (id) {
    case 0:
      return CaseEmpty();
    case 1:
      return CasePoint();
    case 2:
      return CaseLine();
    case 3:
      return CaseSquare();
    case 4:
      return CaseTwoComponents();
    case 5:
      return CaseHole();
    case 6:
      return CaseTouchBorder();
    default:
      return MakeEmpty(1, 1);
  }
}

bool IsAllZero(const BinaryImage &img) {
  return std::ranges::all_of(img.data, [](std::uint8_t v) { return v == 0; });
}

}  // namespace

template <typename TaskT>
class PeryashkinVRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
    if (IsAllZero(input_data_)) {
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

using PeryashkinVRunFuncTestsSEQ = PeryashkinVRunFuncTests<PeryashkinVBinaryComponentContourProcessingSEQ>;
using PeryashkinVRunFuncTestsOMP = PeryashkinVRunFuncTests<PeryashkinVBinaryComponentContourProcessingOMP>;
using PeryashkinVRunFuncTestsTBB = PeryashkinVRunFuncTests<PeryashkinVBinaryComponentContourProcessingTBB>;
using PeryashkinVRunFuncTestsSTL = PeryashkinVRunFuncTests<PeryashkinVBinaryComponentContourProcessingSTL>;

TEST_P(PeryashkinVRunFuncTestsSEQ, BinaryComponentContourSEQ) {
  ExecuteTest(GetParam());
}

TEST_P(PeryashkinVRunFuncTestsOMP, BinaryComponentContourOMP) {
  ExecuteTest(GetParam());
}

TEST_P(PeryashkinVRunFuncTestsTBB, BinaryComponentContourTBB) {
  ExecuteTest(GetParam());
}

TEST_P(PeryashkinVRunFuncTestsSTL, BinaryComponentContourSTL) {
  ExecuteTest(GetParam());
}

namespace {

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(0, "empty"),        std::make_tuple(1, "point"),          std::make_tuple(2, "line"),
    std::make_tuple(3, "square"),       std::make_tuple(4, "two_components"), std::make_tuple(5, "hole"),
    std::make_tuple(6, "touch_border"),
};

const auto kSeqTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<PeryashkinVBinaryComponentContourProcessingSEQ, InType>(
        kTestParam, PPC_SETTINGS_peryashkin_v_binary_component_contour_processing));

const auto kOmpTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<PeryashkinVBinaryComponentContourProcessingOMP, InType>(
        kTestParam, PPC_SETTINGS_peryashkin_v_binary_component_contour_processing));

const auto kTbbTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<PeryashkinVBinaryComponentContourProcessingTBB, InType>(
        kTestParam, PPC_SETTINGS_peryashkin_v_binary_component_contour_processing));

const auto kStlTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<PeryashkinVBinaryComponentContourProcessingSTL, InType>(
        kTestParam, PPC_SETTINGS_peryashkin_v_binary_component_contour_processing));

const auto kSeqValues = ppc::util::ExpandToValues(kSeqTasksList);
const auto kOmpValues = ppc::util::ExpandToValues(kOmpTasksList);
const auto kTbbValues = ppc::util::ExpandToValues(kTbbTasksList);
const auto kStlValues = ppc::util::ExpandToValues(kStlTasksList);

const auto kNameFnSeq = PeryashkinVRunFuncTestsSEQ::PrintFuncTestName<PeryashkinVRunFuncTestsSEQ>;
const auto kNameFnOmp = PeryashkinVRunFuncTestsOMP::PrintFuncTestName<PeryashkinVRunFuncTestsOMP>;
const auto kNameFnTbb = PeryashkinVRunFuncTestsTBB::PrintFuncTestName<PeryashkinVRunFuncTestsTBB>;
const auto kNameFnStl = PeryashkinVRunFuncTestsSTL::PrintFuncTestName<PeryashkinVRunFuncTestsSTL>;

INSTANTIATE_TEST_SUITE_P(FuncTests, PeryashkinVRunFuncTestsSEQ, kSeqValues, kNameFnSeq);
INSTANTIATE_TEST_SUITE_P(FuncTests, PeryashkinVRunFuncTestsOMP, kOmpValues, kNameFnOmp);
INSTANTIATE_TEST_SUITE_P(FuncTests, PeryashkinVRunFuncTestsTBB, kTbbValues, kNameFnTbb);
INSTANTIATE_TEST_SUITE_P(FuncTests, PeryashkinVRunFuncTestsSTL, kStlValues, kNameFnStl);

}  // namespace

}  // namespace peryashkin_v_binary_component_contour_processing
