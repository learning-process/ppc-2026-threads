#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>


#include "krykov_e_sobel_op/common/include/common.hpp"
#include "krykov_e_sobel_op/seq/include/ops_seq.hpp"

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace krykov_e_sobel_op {

class KrykovERunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" +
           std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    int size = std::get<0>(
        std::get<static_cast<size_t>(
            ppc::util::GTestParamIndex::kTestParams)>(GetParam()));

    Image img;
    img.width = size;
    img.height = size;
    img.data.resize(size * size);

    // Генерация тестового изображения:
    // горизонтальный градиент яркости
    for (int y = 0; y < size; ++y) {
      for (int x = 0; x < size; ++x) {
        uint8_t value =
            static_cast<uint8_t>((255 * x) / (size - 1));
        img.data[y * size + x] = {value, value, value};
      }
    }

    input_data_ = img;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (static_cast<int>(output_data.size()) !=
        input_data_.width * input_data_.height) {
      return false;
    }

    int sum = std::accumulate(output_data.begin(),
                              output_data.end(), 0);

    return sum > 0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(KrykovERunFuncTestsThreads, SobelOp) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {
    std::make_tuple(64, "64"),
    std::make_tuple(128, "128"),
    std::make_tuple(256, "256")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<KrykovESobelOpSEQ, InType>(kTestParam, PPC_SETTINGS_krykov_e_sobel_op));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KrykovERunFuncTestsThreads::PrintFuncTestName<KrykovERunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SobelTests, KrykovERunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace krykov_e_sobel_op
