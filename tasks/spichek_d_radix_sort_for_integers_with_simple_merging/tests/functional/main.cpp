#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class RadixSortRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    // Получаем параметры теста: (размер вектора, описание)
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int vector_size = std::get<0>(params);

    // Генерируем случайные целые числа для сортировки
    input_data_.resize(vector_size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-10000, 10000);

    for (int i = 0; i < vector_size; ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Проверяем, что результат действительно отсортирован
    return std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(RadixSortRunFuncTests, RadixSortSequential) {
  ExecuteTest(GetParam());
}

// Определяем наборы данных: (размер массива, описание)
const std::array<TestType, 3> kTestParam = {std::make_tuple(100, "small_vector"),
                                            std::make_tuple(1000, "medium_vector"),
                                            std::make_tuple(5000, "large_vector")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<RadixSortSEQ, InType>(
    kTestParam, PPC_SETTINGS_spichek_d_radix_sort_for_integers_with_simple_merging));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = RadixSortRunFuncTests::PrintFuncTestName<RadixSortRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(RadixSortTests, RadixSortRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
