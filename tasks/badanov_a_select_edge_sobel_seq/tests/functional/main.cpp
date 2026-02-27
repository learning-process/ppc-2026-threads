#include <gtest/gtest.h>

#include <array>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include "badanov_a_select_edge_sobel_seq/common/include/common.hpp"
#include "badanov_a_select_edge_sobel_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace badanov_a_select_edge_sobel_seq {

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
        ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_badanov_a_select_edge_sobel_seq), filename);

    std::ifstream file(abs_path);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + abs_path);
    }

    int width = 0;
    int height = 0;
    file >> width >> height;

    const size_t total_pixels = static_cast<size_t>(width) * static_cast<size_t>(height);
    input_data_.resize(total_pixels);

    for (size_t i = 0; i < total_pixels; ++i) {
      int pixel_value = 0;
      file >> pixel_value;
      input_data_[i] = static_cast<uint8_t>(pixel_value);
    }

    file.close();
  }

  bool CheckAllPixelsZero(const std::vector<uint8_t> &data) const {
    for (const uint8_t pixel : data) {
      if (pixel != 0) {
        return false;
      }
    }
    return true;
  }

  bool CheckImageBorders(const std::vector<uint8_t> &data, int image_width, int image_height) const {
    for (int column = 0; column < image_width; ++column) {
      const size_t index = static_cast<size_t>(column);
      if (data[index] != 0) {
        return false;
      }
    }

    for (int column = 0; column < image_width; ++column) {
      const size_t index = static_cast<size_t>((image_height - 1) * image_width + column);
      if (data[index] != 0) {
        return false;
      }
    }

    for (int row = 0; row < image_height; ++row) {
      const size_t index = static_cast<size_t>(row) * static_cast<size_t>(image_width);
      if (data[index] != 0) {
        return false;
      }
    }

    for (int row = 0; row < image_height; ++row) {
      const size_t index =
          static_cast<size_t>(row) * static_cast<size_t>(image_width) + static_cast<size_t>(image_width - 1);
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

    const int image_width = static_cast<int>(std::sqrt(static_cast<double>(input_data_.size())));
    const int image_height = image_width;

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
        const size_t index = static_cast<size_t>(row) * static_cast<size_t>(image_width) + static_cast<size_t>(column);
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
};

namespace {

TEST_P(BadanovASelectEdgeSobelFuncTests, SobelOnFiles) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(50, "test_1.txt"),  // Простой квадрат
    std::make_tuple(30, "test_2.txt"),  // Градиент
    std::make_tuple(40, "test_3.txt"),  // Диагональная линия
    std::make_tuple(50, "test_4.txt"),  // Пустое изображение
    std::make_tuple(50, "test_6.txt")   // Крест
};

const auto kTestTasksList = ppc::util::AddFuncTask<BadanovASelectEdgeSobelSEQ, InType>(
    kTestParam, PPC_SETTINGS_badanov_a_select_edge_sobel_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BadanovASelectEdgeSobelFuncTests::PrintFuncTestName<BadanovASelectEdgeSobelFuncTests>;

INSTANTIATE_TEST_SUITE_P(SobelEdgeTests, BadanovASelectEdgeSobelFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace badanov_a_select_edge_sobel_seq
