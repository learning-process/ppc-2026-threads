#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "badanov_a_select_edge_sobel/common/include/common.hpp"
#include "badanov_a_select_edge_sobel/omp/include/ops_omp.hpp"
#include "badanov_a_select_edge_sobel/seq/include/ops_seq.hpp"
#include "badanov_a_select_edge_sobel/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace badanov_a_select_edge_sobel {

class BadanovASelectEdgeSobelFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    std::string name = std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
    for (char &c : name) {
      if ((std::isalnum(static_cast<unsigned char>(c)) == 0) && c != '_') {
        c = '_';
      }
    }
    return name;
  }

 protected:
  void SetUp() override {
    const auto &[threshold, filename] =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    threshold_ = threshold;
    filename_ = filename;

    const std::string abs_path =
        ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_badanov_a_select_edge_sobel), filename);

    std::ifstream file(abs_path);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + abs_path);
    }

    int width = 0;
    int height = 0;
    file >> width >> height;

    image_width_ = width;
    image_height_ = height;

    const size_t total_pixels = static_cast<size_t>(width) * static_cast<size_t>(height);
    input_data_.resize(total_pixels);

    for (size_t i = 0; i < total_pixels; ++i) {
      int pixel_value = 0;
      file >> pixel_value;
      input_data_[i] = static_cast<uint8_t>(pixel_value);
    }

    file.close();
  }

  static bool CheckAllPixelsZero(const std::vector<uint8_t> &data) {
    return std::ranges::all_of(data.begin(), data.end(), [](uint8_t pixel) { return pixel == 0; });
  }

  static bool CheckImageBorders(const std::vector<uint8_t> &data, int image_width, int image_height) {
    for (int column = 0; column < image_width; ++column) {
      const auto index = static_cast<size_t>(column);
      if (data[index] != 0) {
        return false;
      }
    }

    for (int column = 0; column < image_width; ++column) {
      const int64_t flat_index =
          (static_cast<int64_t>(image_height - 1) * static_cast<int64_t>(image_width)) + static_cast<int64_t>(column);
      const auto index = static_cast<size_t>(flat_index);
      if (data[index] != 0) {
        return false;
      }
    }

    for (int row = 0; row < image_height; ++row) {
      const auto index = static_cast<size_t>(row) * static_cast<size_t>(image_width);
      if (data[index] != 0) {
        return false;
      }
    }

    for (int row = 0; row < image_height; ++row) {
      const auto index =
          (static_cast<size_t>(row) * static_cast<size_t>(image_width)) + static_cast<size_t>(image_width - 1);
      if (data[index] != 0) {
        return false;
      }
    }

    return true;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.empty()) {
      return false;
    }

    if (output_data.size() != input_data_.size()) {
      return false;
    }

    if (image_height_ < 3 || image_width_ < 3) {
      return output_data == input_data_;
    }

    if (image_width_ != image_height_) {
      return true;
    }

    const int image_width = image_width_;
    const int image_height = image_height_;

    if (!CheckImageBorders(output_data, image_width, image_height)) {
      return false;
    }

    const bool input_all_zeros = CheckAllPixelsZero(input_data_);

    if (input_all_zeros) {
      return CheckAllPixelsZero(output_data);
    }

    bool has_edges = false;
    for (int row = 1; row < image_height - 1 && !has_edges; ++row) {
      for (int column = 1; column < image_width - 1 && !has_edges; ++column) {
        const size_t index =
            (static_cast<size_t>(row) * static_cast<size_t>(image_width)) + static_cast<size_t>(column);
        if (output_data[index] > 0) {
          has_edges = true;
        }
      }
    }

    return has_edges;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  int threshold_{50};
  std::string filename_;
  int image_width_{0};
  int image_height_{0};
};

class BadanovASelectEdgeSobelGradientTests : public ::testing::Test {
 protected:
  void SetUp() override {
    width = 5;
    height = 5;
    input.resize(25);

    for (int row = 0; row < height; ++row) {
      for (int col = 0; col < width; ++col) {
        size_t idx = (row * width) + col;
        input[idx] = static_cast<uint8_t>(row * 50);
      }
    }
  }

  int width{0};
  int height{0};
  InType input;
};

TEST_F(BadanovASelectEdgeSobelGradientTests, GradientComputationProducesEdges) {
  BadanovASelectEdgeSobelOMP task(input);

  task.Validation();
  task.PreProcessing();
  task.Run();
  task.PostProcessing();

  auto output = task.GetOutput();
  ASSERT_FALSE(output.empty());

  bool has_edges = false;
  for (uint8_t pixel : output) {
    if (pixel > 0) {
      has_edges = true;
      break;
    }
  }
  EXPECT_TRUE(has_edges);
}

TEST_F(BadanovASelectEdgeSobelGradientTests, EdgeDetectionWorksWithSharpTransition) {
  InType edge_input(9, 0);
  edge_input[0] = 0;
  edge_input[1] = 0;
  edge_input[2] = 255;
  edge_input[3] = 0;
  edge_input[4] = 0;
  edge_input[5] = 255;
  edge_input[6] = 0;
  edge_input[7] = 0;
  edge_input[8] = 255;

  BadanovASelectEdgeSobelOMP task(edge_input);

  task.Validation();
  task.PreProcessing();
  task.Run();
  task.PostProcessing();

  auto output = task.GetOutput();
  ASSERT_FALSE(output.empty());
  ASSERT_GT(output.size(), static_cast<size_t>(4));
  EXPECT_GT(static_cast<int>(output[4]), 0);
}

// Тесты для граничных случаев
TEST(BadanovASelectEdgeSobelOMPEdgeCases, AllZeroImageOutputAllZero) {
  InType zero_input(100, 0);
  BadanovASelectEdgeSobelOMP task(zero_input);

  task.Validation();
  task.PreProcessing();
  task.Run();
  task.PostProcessing();

  auto output = task.GetOutput();
  for (uint8_t pixel : output) {
    EXPECT_EQ(pixel, 0);
  }
}

TEST(BadanovASelectEdgeSobelOMPEdgeCases, ConstantImageNoEdges) {
  InType constant_input(100, 128);
  BadanovASelectEdgeSobelOMP task(constant_input);

  task.Validation();
  task.PreProcessing();
  task.Run();
  task.PostProcessing();

  auto output = task.GetOutput();
  for (uint8_t pixel : output) {
    EXPECT_EQ(pixel, 0);
  }
}

TEST(BadanovASelectEdgeSobelTBBEdgeCases, VerySmallImage_2x2) {
  InType small_input(4, 0);
  small_input[0] = 0;
  small_input[1] = 255;
  small_input[2] = 255;
  small_input[3] = 0;

  BadanovASelectEdgeSobelTBB task(small_input);
  task.Validation();
  task.PreProcessing();
  task.Run();
  task.PostProcessing();

  auto output = task.GetOutput();
  EXPECT_EQ(output.size(), static_cast<size_t>(4));
  EXPECT_EQ(output, small_input);
}

namespace {

TEST_P(BadanovASelectEdgeSobelFuncTests, SobelOnFiles) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {
    std::make_tuple(50, "test_1.txt"),   // Простой квадрат
    std::make_tuple(30, "test_2.txt"),   // Градиент
    std::make_tuple(40, "test_3.txt"),   // Диагональная линия
    std::make_tuple(50, "test_4.txt"),   // Пустое изображение
    std::make_tuple(50, "test_6.txt"),   // Крест
    std::make_tuple(50, "test_7.txt"),   // Неквадратное 1x8
    std::make_tuple(50, "test_8.txt"),   // Неквадратное 8x1
    std::make_tuple(50, "test_9.txt"),   // Маленькое 2x2
    std::make_tuple(50, "test_10.txt"),  // Маленькое 1x1
    std::make_tuple(50, "test_11.txt")   // 3x2
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BadanovASelectEdgeSobelSEQ, InType>(kTestParam, PPC_SETTINGS_badanov_a_select_edge_sobel),
    ppc::util::AddFuncTask<BadanovASelectEdgeSobelOMP, InType>(kTestParam, PPC_SETTINGS_badanov_a_select_edge_sobel),
    ppc::util::AddFuncTask<BadanovASelectEdgeSobelTBB, InType>(kTestParam, PPC_SETTINGS_badanov_a_select_edge_sobel));
const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BadanovASelectEdgeSobelFuncTests::PrintFuncTestName<BadanovASelectEdgeSobelFuncTests>;

INSTANTIATE_TEST_SUITE_P(SobelEdgeTests, BadanovASelectEdgeSobelFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace badanov_a_select_edge_sobel
