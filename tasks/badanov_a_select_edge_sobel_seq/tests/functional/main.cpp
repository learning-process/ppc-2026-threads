#include <gtest/gtest.h>

#include <array>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include "badanov_a_select_edge_sobel/common/include/common.hpp"
#include "badanov_a_select_edge_sobel/seq/include/ops_seq.hpp"
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
    
    const std::string abs_path =
        ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_badanov_a_select_edge_sobel), "data/" + filename);
    
    std::ifstream file(abs_path);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + abs_path);
    }
    
    int width, height;
    file >> width >> height;
    
    input_data_.resize(width * height);
    for (int i = 0; i < width * height; ++i) {
      int pixel;
      file >> pixel;
      input_data_[i] = static_cast<uint8_t>(pixel);
    }
    
    file.close();
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.empty()) {
      return false;
    }
    
    if (output_data.size() != input_data_.size()) {
      return false;
    }
    
    int width = static_cast<int>(std::sqrt(input_data_.size()));
    int height = width;
    
    for (int col = 0; col < width; ++col) {
      if (output_data[col] != 0) return false;
      if (output_data[(height - 1) * width + col] != 0) return false;
    }
    
    for (int row = 0; row < height; ++row) {
      if (output_data[row * width] != 0) return false;
      if (output_data[row * width + (width - 1)] != 0) return false;
    }
    
    bool all_zeros = true;
    for (uint8_t pixel : input_data_) {
      if (pixel != 0) {
        all_zeros = false;
        break;
      }
    }
    
    if (all_zeros) {
      for (uint8_t pixel : output_data) {
        if (pixel != 0) {
          return false;
        }
      }
      return true;
    }
    
    bool has_edges = false;
    for (int row = 1; row < height - 1 && !has_edges; ++row) {
      for (int col = 1; col < width - 1 && !has_edges; ++col) {
        if (output_data[row * width + col] > 0) {
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
};

namespace {

TEST_P(BadanovASelectEdgeSobelFuncTests, SobelOnFiles) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParam = {
    std::make_tuple(50, "test_1.txt"),   // Простой квадрат
    std::make_tuple(30, "test_2.txt"),   // Градиент
    std::make_tuple(40, "test_3.txt"),   // Диагональная линия
    std::make_tuple(50, "test_4.txt"),   // Пустое изображение
    std::make_tuple(60, "test_5.txt"),   // Шахматная доска
    std::make_tuple(50, "test_6.txt")    // Крест
};

const auto kTestTasksList = ppc::util::AddFuncTask<BadanovASelectEdgeSobel, InType>(
    kTestParam, PPC_SETTINGS_badanov_a_select_edge_sobel);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    BadanovASelectEdgeSobelFuncTests::PrintFuncTestName<BadanovASelectEdgeSobelFuncTests>;

INSTANTIATE_TEST_SUITE_P(SobelEdgeTests, BadanovASelectEdgeSobelFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace badanov_a_select_edge_sobel