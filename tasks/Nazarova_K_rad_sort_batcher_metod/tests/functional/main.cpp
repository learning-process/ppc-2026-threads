#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>

#include "Nazarova_K_rad_sort_batcher_metod/common/include/common.hpp"
#include "Nazarova_K_rad_sort_batcher_metod/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nazarova_k_rad_sort_batcher_metod_processes {

class NazarovaKRadSortBatcherMetodRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType param = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(param);
    expected_ = std::get<1>(param);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_.size()) {
      return false;
    }
    for (std::size_t i = 0; i < output_data.size(); ++i) {
      if (output_data[i] != expected_[i]) {
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
  OutType expected_;
};

namespace {

TEST_P(NazarovaKRadSortBatcherMetodRunFuncTests, RadixSortDoubleBatcherMerge) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {
    TestType{InType{8.8}, OutType{8.8}, "Single"},
    TestType{InType{}, OutType{}, "Empty"},
    TestType{InType{6.6, 3.3}, OutType{3.3, 6.6}, "ReverseTwo"},
    TestType{InType{-0.2, -150.0, -60.5, -4.4, -9.9}, OutType{-150.0, -60.5, -9.9, -4.4, -0.2}, "Negative"},
    TestType{InType{15.7, 0.6, 98.2, 3.75, 7.83, 46.0}, OutType{0.6, 3.75, 7.83, 15.7, 46.0, 98.2}, "Positive"},
    TestType{InType{9e12, 1e3, 5e9, 7e15, 2e11}, OutType{1e3, 5e9, 2e11, 9e12, 7e15}, "Large"},
    TestType{InType{1e-20, 5e-6, 3e-12, 2e-3, 4e-9}, OutType{1e-20, 3e-12, 4e-9, 5e-6, 2e-3}, "Small"},
    TestType{InType{-8.0, -2.0, 0.5, 8.0, 9.0}, OutType{-8.0, -2.0, 0.5, 8.0, 9.0}, "Sorted"},
    TestType{InType{-3.3, 6.6, -10.9, 0.0, 2.2, -1.1}, OutType{-10.9, -3.3, -1.1, 0.0, 2.2, 6.6}, "DifferentSigns"},
    TestType{InType{7.7, 3.3, 7.7, 3.3, 7.7}, OutType{3.3, 3.3, 7.7, 7.7, 7.7}, "Duplicates"},
    TestType{InType{36.6, 25.5, 10.0, 8.9, 6.7, 4.5, 2.2},
             OutType{2.2, 4.5, 6.7, 8.9, 10.0, 25.5, 36.6},
             "ReverseSeven"},
    TestType{InType{0.0, -0.0}, OutType{-0.0, 0.0}, "ZeroAndNegativeZero"}};

const auto kTestTasksList = ppc::util::AddFuncTask<NazarovaKRadSortBatcherMetodSEQ, InType>(
    kTestParam, PPC_SETTINGS_Nazarova_K_rad_sort_batcher_metod);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    NazarovaKRadSortBatcherMetodRunFuncTests::PrintFuncTestName<NazarovaKRadSortBatcherMetodRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(RadixSortBatcherTests, NazarovaKRadSortBatcherMetodRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nazarova_k_rad_sort_batcher_metod_processes
