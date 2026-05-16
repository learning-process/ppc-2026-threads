#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "boltenkov_s_gaussian_kernel/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace boltenkov_s_gaussian_kernel {

class BoltenkovSRunFuncTestsProcesses : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string file_name = params + ".bin";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_boltenkov_s_gaussian_kernel, file_name);
    ReadData(abs_path, input_data_);
    file_name = params + "_output.bin";
    abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_boltenkov_s_gaussian_kernel, file_name);
    ReadData(abs_path, output_data_test_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != std::get<0>(output_data_test_)) {
      return false;
    }
    std::size_t n = std::get<0>(output_data_test_);
    for (std::size_t i = 0; i < n; i++) {
      if (output_data[i].size() != std::get<1>(output_data_test_)) {
        return false;
      }
    }
    std::size_t m = std::get<1>(output_data_test_);
    for (std::size_t i = 0; i < n; i++) {
      for (std::size_t j = 0; j < m; j++) {
        if (output_data[i][j] != std::get<2>(output_data_test_)[i][j]) {
          return false;
        }
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  InType output_data_test_;
  static void ReadData(std::string &abs_path, InType &data) {
    std::ifstream file_stream(abs_path, std::ios::in | std::ios::binary);
    if (!file_stream.is_open()) {
      throw std::runtime_error("Error opening file!\n");
    }
    constexpr std::size_t kMaxSize = 1000;
    int m = -1;
    int n = -1;
    file_stream.read(reinterpret_cast<char *>(&m), sizeof(int));
    file_stream.read(reinterpret_cast<char *>(&n), sizeof(int));
    if (file_stream.fail() || m <= 0 || n <= 0 || std::cmp_greater(n, kMaxSize) || std::cmp_greater(m, kMaxSize)) {
      throw std::runtime_error("invalid input data!\n");
    }
    std::get<0>(data) = static_cast<std::size_t>(n);
    std::get<1>(data) = static_cast<std::size_t>(m);
    std::vector<std::vector<int>> &mtr = std::get<2>(data);
    mtr.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      mtr[i].resize(static_cast<std::size_t>(m));
      file_stream.read(reinterpret_cast<char *>(mtr[i].data()), static_cast<std::streamsize>(sizeof(int) * m));
      if (file_stream.fail()) {
        throw std::runtime_error("Failed to read row " + std::to_string(i));
      }
    }
    file_stream.close();
  }
};

namespace {

TEST_P(BoltenkovSRunFuncTestsProcesses, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 2> kTestParam = {"pic1", "pic2"};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BoltenkovSGaussianKernelSEQ, InType>(kTestParam, PPC_SETTINGS_boltenkov_s_gaussian_kernel));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BoltenkovSRunFuncTestsProcesses::PrintFuncTestName<BoltenkovSRunFuncTestsProcesses>;

INSTANTIATE_TEST_SUITE_P(GaussianKernelTests, BoltenkovSRunFuncTestsProcesses, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace boltenkov_s_gaussian_kernel
