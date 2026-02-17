// tasks/peryashkin_v_binary_component_contour_processing/tests/functional/main.cpp
#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"
#include "peryashkin_v_binary_component_contour_processing/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace peryashkin_v_binary_component_contour_processing {

// TestType в PPC обычно tuple<something, string>
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
  img.data[static_cast<std::size_t>(y) * static_cast<std::size_t>(img.width) + static_cast<std::size_t>(x)] = v;
}

// Набор детерминированных кейсов по id
BinaryImage BuildCase(int id) {
  switch (id) {
    case 0: {  // пустое
      return MakeEmpty(5, 4);
    }
    case 1: {  // точка
      auto im = MakeEmpty(5, 5);
      Set(im, 2, 2, 1);
      return im;
    }
    case 2: {  // линия
      auto im = MakeEmpty(7, 5);
      for (int x = 1; x <= 5; ++x) {
        Set(im, x, 2, 1);
      }
      return im;
    }
    case 3: {  // квадрат 3x3
      auto im = MakeEmpty(6, 6);
      for (int y = 2; y <= 4; ++y) {
        for (int x = 2; x <= 4; ++x) {
          Set(im, x, y, 1);
        }
      }
      return im;
    }
    case 4: {  // 2 компоненты
      auto im = MakeEmpty(8, 6);
      Set(im, 1, 1);
      Set(im, 2, 1);
      Set(im, 1, 2);
      Set(im, 2, 2);
      for (int y = 3; y <= 5; ++y) {
        Set(im, 6, y);
        Set(im, 7, y);
      }
      return im;
    }
    case 5: {  // “дырка” (кольцо 5x5)
      auto im = MakeEmpty(7, 7);
      for (int x = 1; x <= 5; ++x) {
        Set(im, x, 1);
        Set(im, x, 5);
      }
      for (int y = 1; y <= 5; ++y) {
        Set(im, 1, y);
        Set(im, 5, y);
      }
      return im;
    }
    case 6: {  // касается границы
      auto im = MakeEmpty(5, 5);
      for (int y = 0; y <= 2; ++y) {
        for (int x = 0; x <= 1; ++x) {
          Set(im, x, y, 1);
        }
      }
      return im;
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
    // Здесь проверяем инварианты, не “рандом”
    // 1) детерминированность: для пустого => 0 компонент
    if (input_data_.data.empty()) {
      return false;
    }
    if (std::all_of(input_data_.data.begin(), input_data_.data.end(), [](uint8_t v) { return v == 0; })) {
      return output_data.empty();
    }
    // 2) для непустых кейсов — хотя бы 1 hull
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

// ---- Доп. unit-тесты для покрытия веток Validation/edge ----

TEST(PeryashkinVBinaryComponentContourProcessingSEQ_Unit, ValidationFailsOnBadSizes) {
  BinaryImage bad;
  bad.width = 4;
  bad.height = 4;
  bad.data.assign(3, 1);  // неверный размер
  PeryashkinVBinaryComponentContourProcessingSEQ task(bad);
  EXPECT_FALSE(task.Validation());
}

TEST(PeryashkinVBinaryComponentContourProcessingSEQ_Unit, ValidationFailsOnNonPositiveDims) {
  BinaryImage bad;
  bad.width = 0;
  bad.height = 5;
  bad.data.clear();
  PeryashkinVBinaryComponentContourProcessingSEQ task(bad);
  EXPECT_FALSE(task.Validation());
}

TEST(PeryashkinVBinaryComponentContourProcessingSEQ_Unit, PipelineReturnsFalseIfInvalid) {
  BinaryImage bad;
  bad.width = 3;
  bad.height = 3;
  bad.data.assign(1, 1);
  PeryashkinVBinaryComponentContourProcessingSEQ task(bad);
  EXPECT_FALSE(task.Validation());
  EXPECT_TRUE(task.PreProcessing());  // PreProcessing не обязан падать
  EXPECT_FALSE(task.Run());           // RunImpl проверяет валидность
}

}  // namespace peryashkin_v_binary_component_contour_processing
