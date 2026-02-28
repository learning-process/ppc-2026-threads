#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "ivanova_p_marking_components_on_binary_image/data/image_generator.hpp"
#include "ivanova_p_marking_components_on_binary_image/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  IvanovaPRunFuncTestsThreads() = default;

  static std::string PrintTestParam(const TestType &test_param) {
    return "test_" + std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

  static int GetExpectedComponents(int test_case) {
    switch (test_case) {
      case 1:
      case 4:
      case 6:
      case 9:
      case 10:
      case 11:
        return 1;
      case 2:
      case 12:
      case 13:
        return 2;
      case 3:
        return 3;
      case 5:
        return 0;
      case 7:
        return 4;
      case 8:
      case 14:
        return 9;
      default:
        return 0;
    }
  }

 protected:
  void SetUp() override {
    // GetParam() возвращает FuncTestParam, который содержит (функция, имя, TestType)
    // TestType это std::tuple<int, std::string>
    const TestType &test_param = std::get<2>(GetParam());
    current_test_case_ = std::get<0>(test_param);

    // Определяем, является ли это файловым тестом (тесты 11-14)
    is_file_test_ = (current_test_case_ >= 11 && current_test_case_ <= 14);
  }

  // main.cpp (часть с CheckTestOutputData)
  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() < 3) {
      return false;
    }

    int out_width = output_data[0];
    int out_height = output_data[1];
    int num_components = output_data[2];

    // Проверяем базовые параметры
    if (out_width != test_image.width || out_height != test_image.height) {
      return false;
    }

    // Проверяем количество компонент
    int expected_components = GetExpectedComponents(current_test_case_);
    if (num_components != expected_components) {
      return false;
    }

    // Проверяем корректность меток
    std::vector<bool> found_labels(static_cast<size_t>(num_components) + 1, false);

    for (size_t i = 3; i < output_data.size(); ++i) {
      int label = output_data[i];
      size_t idx = i - 3;

      if (idx >= test_image.data.size()) {
        return false;
      }

      uint8_t original_pixel = test_image.data[idx];

      if (!ValidatePixel(label, original_pixel, num_components, found_labels)) {
        return false;
      }
    }

    return AreAllLabelsUsed(found_labels, num_components);
  }

  InType GetTestInputData() final {
    // Инициализируем изображение здесь, перед созданием задачи
    if (is_file_test_) {
      // Загружаем изображение из файла
      std::string filename;
      switch (current_test_case_) {
        case 11:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image.txt";
          break;
        case 12:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image2.txt";
          break;
        case 13:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image3.txt";
          break;
        case 14:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image4.txt";
          break;
        default:
          filename = "";
      }
      test_image = LoadImageFromTxt(filename);
    } else {
      // Создаем тестовое изображение программно
      const int width = 100;
      const int height = 100;
      test_image = CreateTestImage(width, height, current_test_case_);
    }

    return current_test_case_;
  }

 private:
  int current_test_case_ = 0;
  bool is_file_test_ = false;

  static bool ValidatePixel(int label, uint8_t original_pixel, int num_components, std::vector<bool> &found_labels) {
    if (original_pixel == 0) {
      return label == 0;
    }

    if (label < 1 || label > num_components) {
      return false;
    }

    found_labels[static_cast<size_t>(label)] = true;
    return true;
  }

  [[nodiscard]] static bool AreAllLabelsUsed(const std::vector<bool> &found_labels, int num_components) {
    for (int i = 1; i <= num_components; ++i) {
      if (!found_labels[static_cast<size_t>(i)]) {
        return false;
      }
    }
    return true;
  }
};

namespace {

TEST_P(IvanovaPRunFuncTestsThreads, MarkingComponentsTest) {
  ExecuteTest(GetParam());
}

// Тестовые случаи: (код_теста, описание)
const std::array<TestType, 14> kTestParam = {
    std::make_tuple(1, "single_component"),   std::make_tuple(2, "two_components"),
    std::make_tuple(3, "three_components"),   std::make_tuple(4, "connected_components"),
    std::make_tuple(5, "all_background"),     std::make_tuple(6, "single_pixel"),
    std::make_tuple(7, "diagonal_neighbors"), std::make_tuple(8, "many_small_components"),
    std::make_tuple(9, "horizontal_line"),    std::make_tuple(10, "vertical_line"),
    std::make_tuple(11, "file_image"),        std::make_tuple(12, "file_image2"),
    std::make_tuple(13, "file_image3"),       std::make_tuple(14, "file_image4")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<IvanovaPMarkingComponentsOnBinaryImageSEQ, InType>(
    kTestParam, PPC_SETTINGS_ivanova_p_marking_components_on_binary_image));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = IvanovaPRunFuncTestsThreads::PrintFuncTestName<IvanovaPRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(ImageTests, IvanovaPRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace ivanova_p_marking_components_on_binary_image
