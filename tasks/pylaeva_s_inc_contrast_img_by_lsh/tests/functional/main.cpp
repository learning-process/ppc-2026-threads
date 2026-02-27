#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"
#include "pylaeva_s_inc_contrast_img_by_lsh/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

class PylaevaSRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param) + "_image_with_size_" + std::to_string(std::get<0>(test_param).size());
  }

 protected:
  void SetUp() override {
    test_params_ = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(test_params_);
    expected_data_ = std::get<1>(test_params_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_data_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  TestType test_params_;
  OutType expected_data_;
};

namespace {

// Тестовые данные
// 1. Простое изображение с градиентом
const std::vector<uint8_t> kImage1 = {0, 50, 100, 150, 200, 255};
const std::vector<uint8_t> kExpected1 = {0, 50, 100, 150, 200, 255};

// 2. Изображение с узким диапазоном значений
const std::vector<uint8_t> kImage2 = {100, 110, 120, 130, 140, 150};
const std::vector<uint8_t> kExpected2 = {0, 51, 102, 153, 204, 255};

// 3. Изображение с очень узким диапазоном
const std::vector<uint8_t> kImage3 = {200, 201, 202, 203, 204, 205};
const std::vector<uint8_t> kExpected3 = {0, 51, 102, 153, 204, 255};

// 4. Изображение с постоянными значениями
const std::vector<uint8_t> kImage4 = {128, 128, 128, 128, 128, 128};
const std::vector<uint8_t> kExpected4 = {128, 128, 128, 128, 128, 128};

// 5. Изображение с экстремальными значениями
const std::vector<uint8_t> kImage5 = {0, 255, 0, 255, 0, 255};
const std::vector<uint8_t> kExpected5 = {0, 255, 0, 255, 0, 255};

// 6. Изображение с неравномерным распределением
const std::vector<uint8_t> kImage6 = {10, 10, 10, 10, 10, 200};
const std::vector<uint8_t> kExpected6 = {0, 0, 0, 0, 0, 255};

// 7. Изображение с двумя пиками
const std::vector<uint8_t> kImage7 = {30, 30, 30, 180, 180, 180};
const std::vector<uint8_t> kExpected7 = {0, 0, 0, 255, 255, 255};

// 8. Изображение с промежуточными значениями и шумом
const std::vector<uint8_t> kImage8 = {45, 67, 89, 123, 156, 178, 200, 210};
const std::vector<uint8_t> kExpected8 = {0, 34, 68, 121, 172, 206, 240, 255};

// 9. Изображение с очень низкими значениями (темное)
const std::vector<uint8_t> kImage9 = {10, 15, 20, 25, 30, 35, 40, 45};
const std::vector<uint8_t> kExpected9 = {0, 36, 73, 109, 146, 182, 219, 255};

// 10. Изображение с очень высокими значениями (светлое)
const std::vector<uint8_t> kImage10 = {210, 215, 220, 225, 230, 235, 240, 245};
const std::vector<uint8_t> kExpected10 = {0, 36, 73, 109, 146, 182, 219, 255};

// 11. КРАЙНИЙ СЛУЧАЙ: Изображение с одним элементом (минимальное)
const std::vector<uint8_t> kImage11 = {100};
const std::vector<uint8_t> kExpected11 = {100};  // Должно остаться без изменений

// 12. КРАЙНИЙ СЛУЧАЙ: Изображение с одним элементом (максимальное)
const std::vector<uint8_t> kImage12 = {255};
const std::vector<uint8_t> kExpected12 = {255};

// 13. КРАЙНИЙ СЛУЧАЙ: Изображение с одним элементом (минимальное)
const std::vector<uint8_t> kImage13 = {0};
const std::vector<uint8_t> kExpected13 = {0};

// 14. КРАЙНИЙ СЛУЧАЙ: Два элемента с минимальной разницей
const std::vector<uint8_t> kImage14 = {100, 101};
const std::vector<uint8_t> kExpected14 = {0, 255};

// 15. КРАЙНИЙ СЛУЧАЙ: Два элемента с максимальной разницей
const std::vector<uint8_t> kImage15 = {0, 255};
const std::vector<uint8_t> kExpected15 = {0, 255};

// 16. КРАЙНИЙ СЛУЧАЙ: Все значения на нижней границе
const std::vector<uint8_t> kImage16 = {0, 0, 0, 0, 0, 0};
const std::vector<uint8_t> kExpected16 = {0, 0, 0, 0, 0, 0};

// 17. КРАЙНИЙ СЛУЧАЙ: Все значения на верхней границе
const std::vector<uint8_t> kImage17 = {255, 255, 255, 255, 255, 255};
const std::vector<uint8_t> kExpected17 = {255, 255, 255, 255, 255, 255};

TEST_P(PylaevaSRunFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

// Создаем массив тестовых параметров
const std::array<TestType, 17> kTestParam = {std::make_tuple(kImage1, kExpected1, "gradient"),
                                             std::make_tuple(kImage2, kExpected2, "narrow_range"),
                                             std::make_tuple(kImage3, kExpected3, "very_narrow_range"),
                                             std::make_tuple(kImage4, kExpected4, "constant"),
                                             std::make_tuple(kImage5, kExpected5, "extreme"),
                                             std::make_tuple(kImage6, kExpected6, "uneven"),
                                             std::make_tuple(kImage7, kExpected7, "two_peaks"),
                                             std::make_tuple(kImage8, kExpected8, "noisy"),
                                             std::make_tuple(kImage9, kExpected9, "dark"),
                                             std::make_tuple(kImage10, kExpected10, "bright"),
                                             std::make_tuple(kImage11, kExpected11, "single_element_min"),
                                             std::make_tuple(kImage12, kExpected12, "single_element_max"),
                                             std::make_tuple(kImage13, kExpected13, "single_element_zero"),
                                             std::make_tuple(kImage14, kExpected14, "two_elements_min_diff"),
                                             std::make_tuple(kImage15, kExpected15, "two_elements_max_diff"),
                                             std::make_tuple(kImage16, kExpected16, "all_zero"),
                                             std::make_tuple(kImage17, kExpected17, "all_max")};
const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<PylaevaSIncContrastImgByLshSEQ, InType>(
    kTestParam, PPC_SETTINGS_pylaeva_s_inc_contrast_img_by_lsh));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = PylaevaSRunFuncTestsThreads::PrintFuncTestName<PylaevaSRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PylaevaSRunFuncTests, PylaevaSRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
