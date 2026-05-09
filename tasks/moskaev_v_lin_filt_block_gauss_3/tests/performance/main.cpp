#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"
#include "moskaev_v_lin_filt_block_gauss_3/omp/include/ops_omp.hpp"
#include "moskaev_v_lin_filt_block_gauss_3/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

class MoskaevVLinFiltBlockGauss3PerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    int width = 2048;
    int height = 2048;
    int channels = 3;
    int block_size = 64;

    std::vector<uint8_t> image_data(static_cast<size_t>(width) * static_cast<size_t>(height) *
                                    static_cast<size_t>(channels));

    for (int row = 0; row < height; ++row) {
      for (int col = 0; col < width; ++col) {
        for (int channel = 0; channel < channels; ++channel) {
          const size_t idx = ((static_cast<size_t>(row) * static_cast<size_t>(width) + static_cast<size_t>(col)) *
                              static_cast<size_t>(channels)) +
                             static_cast<size_t>(channel);
          image_data[idx] = static_cast<uint8_t>((row + col + channel * 85) % 256);
        }
      }
    }

    input_data_ = std::make_tuple(width, height, channels, block_size, image_data);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == std::get<4>(input_data_).size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(MoskaevVLinFiltBlockGauss3PerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = std::tuple_cat(
    ppc::util::MakeAllPerfTasks<InType, MoskaevVLinFiltBlockGauss3SEQ>(PPC_SETTINGS_moskaev_v_lin_filt_block_gauss_3),
    ppc::util::MakeAllPerfTasks<InType, MoskaevVLinFiltBlockGauss3OMP>(PPC_SETTINGS_moskaev_v_lin_filt_block_gauss_3));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MoskaevVLinFiltBlockGauss3PerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTests, MoskaevVLinFiltBlockGauss3PerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace moskaev_v_lin_filt_block_gauss_3
