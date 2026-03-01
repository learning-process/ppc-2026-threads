#include <gtest/gtest.h>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"
#include "chaschin_vladimir_linear_image_filtration_seq/mpi/include/ops_mpi.hpp"
#include "chaschin_vladimir_linear_image_filtration_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace chaschin_v_linear_image_filtration_seq {

class ChaschinVRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kCount = 512;

  void SetUp() override {
    const int width = kCount;
    const int height = kCount;

    // ---------- deterministic image ----------
    std::vector<float> image(width * height);
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        image[i * width + j] = static_cast<float>((i * 131 + j * 17) % 256);
      }
    }

    input_data_ = std::make_tuple(image, width, height);

    // ---------- Gaussian kernel ----------
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
          const int ni = i + di;
          if (ni < 0 || ni >= height) {
            continue;
          }

          for (int dj = -1; dj <= 1; ++dj) {
            const int nj = j + dj;
            if (nj < 0 || nj >= width) {
              continue;
            }

            acc += image[ni * width + nj] * k[di + 1][dj + 1];
          }
        }

        expected_output_[i * width + j] = acc;
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }

    constexpr float eps = 1e-4f;
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

TEST_P(ChaschinVRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, NesterovATestTaskMPI, ChaschinVLinearFiltrationSEQ>(
    PPC_SETTINGS_example_processes);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ChaschinVRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ChaschinVRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chaschin_v_linear_image_filtration_seq
