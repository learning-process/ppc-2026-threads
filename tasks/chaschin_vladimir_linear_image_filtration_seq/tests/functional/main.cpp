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

#include "example_processes/common/include/common.hpp"
#include "example_processes/mpi/include/ops_mpi.hpp"
#include "example_processes/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace chaschin_v_linear_image_filtration_seq {

class ChaschinVRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType& test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int size = std::get<0>(params);

    const int width = size;
    const int height = size;

    // ---------- deterministic image ----------
    std::vector<float> image(width * height);
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        image[i * width + j] = static_cast<float>((i + 1) * (j + 3) % 256);
      }
    }

    input_data_ = std::make_tuple(image, width, height);

    // ---------- Gaussian kernel 3x3 ----------
    const float k[3][3] = {
        {1.f / 16.f, 2.f / 16.f, 1.f / 16.f},
        {2.f / 16.f, 4.f / 16.f, 2.f / 16.f},
        {1.f / 16.f, 2.f / 16.f, 1.f / 16.f},
    };

    expected_output_.resize(width * height, 0.0f);

    // ---------- reference convolution ----------
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        float acc = 0.0f;

        for (int di = -1; di <= 1; ++di) {
          for (int dj = -1; dj <= 1; ++dj) {
            const int ni = i + di;
            const int nj = j + dj;

            if (ni >= 0 && ni < height && nj >= 0 && nj < width) {
              acc += image[ni * width + nj] * k[di + 1][dj + 1];
            }
          }
        }

        expected_output_[i * width + j] = acc;
      }
    }
  }

   bool CheckTestOutputData(OutType& output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }

    constexpr float eps = 1e-5f;
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::fabs(output_data[i] - expected_output_[i]) > eps) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
    OutType expected_output_;
};

namespace {

TEST_P(ChaschinVRunFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(4, "4"),
    std::make_tuple(8, "8"),
    std::make_tuple(16, "16"),
    std::make_tuple(32, "32"),
    std::make_tuple(64, "64"),
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<NesterovATestTaskMPI, InType>(kTestParam, PPC_SETTINGS_example_processes),
    ppc::util::AddFuncTask<ChaschinVLinearFiltrationSEQ, InType>(kTestParam, PPC_SETTINGS_example_processes));

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<ChaschinVLinearFiltrationSEQ, InType>(
        kTestParam, PPC_SETTINGS_example_processes));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    ChaschinVRunFuncTests::PrintFuncTestName<ChaschinVRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(
    LinearGaussianTests,
    ChaschinVRunFuncTests,
    kGtestValues,
    kPerfTestName);

}  // namespace

}  // namespace chaschin_v_linear_image_filtration_seq
