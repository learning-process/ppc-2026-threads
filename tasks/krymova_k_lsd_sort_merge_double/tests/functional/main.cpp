#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>

#include "krymova_k_lsd_sort_merge_double/common/include/common.hpp"
#include "krymova_k_lsd_sort_merge_double/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"  // для ppc::util::GTestParamIndex

namespace krymova_k_lsd_sort_merge_double {

class KrymovaKFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);
    std::string type = std::get<1>(params);

    std::random_device rd;
    std::mt19937 gen(rd());

    input_data_.resize(size);

    if (type == "empty") {
      input_data_.clear();
    } else if (type == "single") {
      input_data_.resize(1);
      input_data_[0] = 42.0;
    } else if (type == "random") {
      std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
      for (int i = 0; i < size; ++i) {
        input_data_[i] = dist(gen);
      }
    } else if (type == "sorted") {
      for (int i = 0; i < size; ++i) {
        input_data_[i] = static_cast<double>(i);
      }
    } else if (type == "reverse") {
      for (int i = 0; i < size; ++i) {
        input_data_[i] = static_cast<double>(size - i);
      }
    } else if (type == "constant") {
      std::ranges::fill(input_data_, 42.0);
    } else if (type == "negative") {
      std::uniform_real_distribution<double> dist(-1000.0, -1.0);
      for (int i = 0; i < size; ++i) {
        input_data_[i] = dist(gen);
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Проверяем, что массив отсортирован
    for (size_t i = 1; i < output_data.size(); ++i) {
      if (output_data[i] < output_data[i - 1]) {
        return false;
      }
    }

    // Проверяем, что размер не изменился
    if (input_data_.size() != output_data.size()) {
      return false;
    }

    // Для непустых массивов проверяем, что все элементы на месте
    if (!input_data_.empty() && !output_data.empty()) {
      OutType input_copy = input_data_;
      const OutType& output_copy = output_data;

      std::ranges::sort(input_copy);

      for (size_t i = 0; i < input_copy.size(); ++i) {
        if (std::abs(input_copy[i] - output_copy[i]) > 1e-10) {
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
};

namespace {

TEST_P(KrymovaKFuncTests, TestSorting) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 14> kTestParam = {
    std::make_tuple(0, "empty"),         std::make_tuple(1, "single"),
    std::make_tuple(10, "random_small"), std::make_tuple(100, "random_medium"),
    std::make_tuple(1000, "random_large"), std::make_tuple(10000, "random_very_large"),
    std::make_tuple(100, "sorted"),       std::make_tuple(100, "reverse"),
    std::make_tuple(100, "constant"),     std::make_tuple(100, "negative"),
    std::make_tuple(1000, "negative_large"), std::make_tuple(10, "mixed"),
    std::make_tuple(100, "mixed"),        std::make_tuple(100000, "random_huge")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<KrymovaKLsdSortMergeDoubleSEQ, InType>(
    kTestParam, PPC_SETTINGS_krymova_k_lsd_sort_merge_double));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KrymovaKFuncTests::PrintFuncTestName<KrymovaKFuncTests>;

INSTANTIATE_TEST_SUITE_P(LsdSortMergeTests, KrymovaKFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace krymova_k_lsd_sort_merge_double