#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "ivanova_p_marking_components_on_binary_image/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  IvanovaPRunFuncTestsThreads() : current_test_case_(0), is_file_test_(false) {}

  static std::string PrintTestParam(const TestType &test_param) {
    return "test_" + std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
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
    int expected_components = 0;

    switch (current_test_case_) {
      case 1:
        expected_components = 1;
        break;
      case 2:
        expected_components = 2;
        break;
      case 3:
        expected_components = 3;
        break;
      case 4:
        expected_components = 1;
        break;
      case 5:
        expected_components = 0;
        break;  // Все фон
      case 6:
        expected_components = 1;
        break;  // Одиночный пиксель
      case 7:
        expected_components = 4;
        break;  // Диагональные соседи - 4 компоненты
      case 8:
        expected_components = 9;
        break;  // 9 маленьких компонент
      case 9:
        expected_components = 1;
        break;  // Горизонтальная линия
      case 10:
        expected_components = 1;
        break;  // Вертикальная линия
      case 11:
        expected_components = 1;
        break;  // image.txt - 1 компонента
      case 12:
        expected_components = 2;
        break;  // image2.txt - 2 компоненты
      case 13:
        expected_components = 2;
        break;  // image3.txt - 2 компоненты
      case 14:
        expected_components = 9;
        break;  // image4.txt - 9 компонент
      default:
        expected_components = num_components;
    }

    if (num_components != expected_components) {
      return false;
    }

    // Проверяем корректность меток
    std::vector<bool> found_labels(num_components + 1, false);

    for (size_t i = 3; i < output_data.size(); ++i) {
      int label = output_data[i];
      int idx = static_cast<int>(i - 3);
      uint8_t original_pixel = test_image.data[idx];

      if (original_pixel == 0) {
        if (label != 0) {
          return false;
        }
      } else {
        if (label < 1 || label > num_components) {
          return false;
        }
        found_labels[label] = true;
      }
    }

    // Проверяем, что все метки использованы
    for (int i = 1; i <= num_components; ++i) {
      if (!found_labels[i]) {
        return false;
      }
    }

    return true;
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
