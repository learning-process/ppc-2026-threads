#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <fstream>
#include <string>

#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"
#include "rozenberg_a_quicksort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace rozenberg_a_quicksort_simple_merge {

class RozenbergARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    input_data_.clear();
    output_data_.clear();
    TestType filename =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()) + ".txt";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_rozenberg_a_quicksort_simple_merge, filename);
    std::ifstream file(abs_path);

    if (file.is_open()) {
      int size = 0;
      file >> size;

      InType input_data(size);
      for (int i = 0; i < size; i++) {
        file >> input_data[i];
      }

      OutType output_data(size);
      for (int i = 0; i < size; i++) {
        file >> output_data[i];
      }
      input_data_ = input_data;
      output_data_ = output_data;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType output_data_;
};

namespace {

TEST_P(RozenbergARunFuncTestsThreads, Quicksort) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 9> kTestParam = {"basic_test",       "duplicate_test",       "inorder_test",
                                            "mixed_sign_test",  "negative_values_test", "random_data_test_2",
                                            "random_data_test", "reverse_order_test",   "single_element_test"};

const auto kTestTasksList = ppc::util::AddFuncTask<RozenbergAQuicksortSimpleMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_rozenberg_a_quicksort_simple_merge);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = RozenbergARunFuncTestsThreads::PrintFuncTestName<RozenbergARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(QuicksortSimpleMergeTests, RozenbergARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace rozenberg_a_quicksort_simple_merge
