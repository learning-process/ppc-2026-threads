#include <gtest/gtest.h>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "boltenkov_s_gaussian_kernel/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace boltenkov_s_gaussian_kernel {

class BoltenkovSRunPerfTestProcesses : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void readData(std::string &abs_path, InType &data) {
    std::ifstream file_stream(abs_path, std::ios::in | std::ios::binary);
    if (!file_stream.is_open()) {
      throw std::runtime_error("Error opening file!\n");
    }
    int m = -1;
    int n = -1;
    file_stream.read(reinterpret_cast<char *>(&m), sizeof(int));
    file_stream.read(reinterpret_cast<char *>(&n), sizeof(int));
    if (m <= 0 || n <= 0) {
      throw std::runtime_error("invalid input data!\n");
    }
    std::get<0>(data) = static_cast<std::size_t>(n);
    std::get<1>(data) = static_cast<std::size_t>(m);
    std::vector<std::vector<int>> &mtr = std::get<2>(data);
    mtr.resize(static_cast<std::size_t>(n), std::vector<int>(static_cast<std::size_t>(m)));
    for (int i = 0; i < n; i++) {
      file_stream.read(reinterpret_cast<char *>(mtr[i].data()), static_cast<std::streamsize>(sizeof(int) * m));
    }
    file_stream.close();
  }

  void SetUp() override {
    std::string file_name = "pic3.bin";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_boltenkov_s_gaussian_kernel, file_name);
    readData(abs_path, input_data_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != std::get<0>(input_data_)) {
      return false;
    }
    std::size_t n = std::get<0>(input_data_);
    for (std::size_t i = 0; i < n; i++) {
      if (output_data[i].size() != std::get<1>(input_data_)) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(BoltenkovSRunPerfTestProcesses, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BoltenkovSGaussianKernelSEQ>(PPC_SETTINGS_boltenkov_s_gaussian_kernel);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BoltenkovSRunPerfTestProcesses::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BoltenkovSRunPerfTestProcesses, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace boltenkov_s_gaussian_kernel
