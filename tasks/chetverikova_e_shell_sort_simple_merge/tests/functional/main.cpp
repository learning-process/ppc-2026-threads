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

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"
#include "chetverikova_e_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

class ChetverikovaERunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
  private:
    InType input_data_;
    OutType expected_data_;
    int n{};
  public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string filename = params + ".txt";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_chetverikova_e_shell_sort_simple_merge, filename);
    std::ifstream file(abs_path);
    if (!file.is_open()) {
      throw std::runtime_error("Failed to open file: " + filename);
    }
    
    if (!(file >> n)) {
      throw std::runtime_error("Failed to read required parameters");
    }
    input_data_.resize(n);
    expected_data_.resize(n);

    for (int i = 0; i < n; ++i) {
        if (!(file >> input_data_[i])) {
            throw std::runtime_error("Failed to read input data at index " + std::to_string(i));
        }
    }

    for (int i = 0; i < n; ++i) {
        if (!(file >> expected_data_[i])) {
            throw std::runtime_error("Failed to read expected data at index " + std::to_string(i));
        }
    }
    file.close();
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

};

namespace {

TEST_P(ChetverikovaERunFuncTestsThreads, ShellSortMergeTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::string("test1"), std::string("test2"), std::string("test3")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<ChetverikovaEShellSortSimpleMergeSEQ, InType>(kTestParam, PPC_SETTINGS_chetverikova_e_shell_sort_simple_merge));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ChetverikovaERunFuncTestsThreads::PrintFuncTestName<ChetverikovaERunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(ShellSortMergeTests, ChetverikovaERunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chetverikova_e_shell_sort_simple_merge
