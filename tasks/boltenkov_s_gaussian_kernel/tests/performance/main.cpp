#include <gtest/gtest.h>

#include <cstddef>
#include <fstream>
#include <ios>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "boltenkov_s_gaussian_kernel/all/include/ops_all.hpp"
#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "boltenkov_s_gaussian_kernel/omp/include/ops_omp.hpp"
#include "boltenkov_s_gaussian_kernel/seq/include/ops_seq.hpp"
#include "boltenkov_s_gaussian_kernel/stl/include/ops_stl.hpp"
#include "boltenkov_s_gaussian_kernel/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace boltenkov_s_gaussian_kernel {

class BoltenkovSRunPerfTestProcesses : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  static void ReadData(std::string &abs_path, InType &data) {
    std::ifstream file_stream(abs_path, std::ios::in | std::ios::binary);
    if (!file_stream.is_open()) {
      throw std::runtime_error("Error opening file!\n");
    }

    const int k_min_size_for_gen = 1000;
    const int k_max_size = 2000;
    int m = -1;
    int n = -1;
    file_stream.read(reinterpret_cast<char *>(&m), sizeof(int));
    file_stream.read(reinterpret_cast<char *>(&n), sizeof(int));
    if (file_stream.fail()) {
      throw std::runtime_error("Failed to read matrix dimensions from file");
    }

    if (m <= 0 || n <= 0 || n > k_max_size || m > k_max_size) {
      throw std::runtime_error("Matrix dimensions exceed maximum allowed size (2000)");
    }

    if (n < k_min_size_for_gen || m < k_min_size_for_gen) {
      n = k_max_size;
      m = k_max_size;
      std::get<0>(data) = static_cast<std::size_t>(n);
      std::get<1>(data) = static_cast<std::size_t>(m);
      std::vector<std::vector<int>> &mtr = std::get<2>(data);
      mtr.resize(static_cast<std::size_t>(n));

      std::random_device rd;
      std::mt19937 gen(rd());

      std::uniform_int_distribution<int> dist(0, 255);

      for (int i = 0; i < n; ++i) {
        mtr[i].resize(static_cast<std::size_t>(m));
        for (int j = 0; j < m; ++j) {
          mtr[i][j] = dist(gen);
        }
      }
      file_stream.close();
      return;
    }

    std::get<0>(data) = static_cast<std::size_t>(n);
    std::get<1>(data) = static_cast<std::size_t>(m);
    std::vector<std::vector<int>> &mtr = std::get<2>(data);
    mtr.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      mtr[i].resize(static_cast<std::size_t>(m));
      file_stream.read(reinterpret_cast<char *>(mtr[i].data()), static_cast<std::streamsize>(sizeof(int) * m));
      if (file_stream.fail()) {
        throw std::runtime_error("Failed to read row " + std::to_string(i));
      }
    }
    file_stream.close();
  }

  void SetUp() override {
    std::string file_name = "pic3.bin";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_boltenkov_s_gaussian_kernel, file_name);
    ReadData(abs_path, input_data_);
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
    ppc::util::MakeAllPerfTasks<InType, BoltenkovSGaussianKernelSEQ, BoltenkovSGaussianKernelOMP,
                                BoltenkovSGaussianKernelTBB, BoltenkovSGaussianKernelSTL, BoltenkovSGaussianKernelALL>(
        PPC_SETTINGS_boltenkov_s_gaussian_kernel);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BoltenkovSRunPerfTestProcesses::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BoltenkovSRunPerfTestProcesses, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace boltenkov_s_gaussian_kernel
