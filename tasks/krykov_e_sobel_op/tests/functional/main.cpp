#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "krykov_e_sobel_op/common/include/common.hpp"
#include "krykov_e_sobel_op/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace krykov_e_sobel_op {

class KrykovERunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    int test_id = std::get<0>(std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));

    const int size = 3;

    Image img;
    img.width = size;
    img.height = size;
    img.data.resize(9);

    expected_output_.assign(9, 0);

    switch (test_id) {
      case 0: {
        // Все пиксели = 100
        for (auto &p : img.data) {
          p = {100, 100, 100};
        }
        // Собель = 0
        break;
      }

      case 1: {
        // Вертикальная граница
        // 0 | 255
        img.data = {{0, 0, 0},       {255, 255, 255}, {255, 255, 255}, {0, 0, 0},      {255, 255, 255},
                    {255, 255, 255}, {0, 0, 0},       {255, 255, 255}, {255, 255, 255}};

        // Центр (1,1)
        // Gx = (-1*0 + 1*255)
        //    + (-2*0 + 2*255)
        //    + (-1*0 + 1*255)
        //    = 255 + 510 + 255 = 1020
        expected_output_[4] = 1020;
        break;
      }

      case 2: {
        // Горизонтальная граница
        img.data = {{0, 0, 0},       {0, 0, 0},       {0, 0, 0},       {255, 255, 255}, {255, 255, 255},
                    {255, 255, 255}, {255, 255, 255}, {255, 255, 255}, {255, 255, 255}};

        // Gy аналогично = 1020
        expected_output_[4] = 1020;
        break;
      }
    }

    input_data_ = img;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_output_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(KrykovERunFuncTestsThreads, SobelOp) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(0, "Constant3x3"), std::make_tuple(1, "VerticalEdge3x3"),
                                            std::make_tuple(2, "HorizontalEdge3x3")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<KrykovESobelOpSEQ, InType>(kTestParam, PPC_SETTINGS_krykov_e_sobel_op));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = KrykovERunFuncTestsThreads::PrintFuncTestName<KrykovERunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SobelTests, KrykovERunFuncTestsThreads, kGtestValues, kTestName);

}  // namespace
}  // namespace krykov_e_sobel_op
