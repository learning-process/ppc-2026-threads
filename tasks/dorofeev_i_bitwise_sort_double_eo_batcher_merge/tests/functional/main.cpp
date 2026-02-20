#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/all/include/ops_all.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/omp/include/ops_omp.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/seq/include/ops_seq.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/stl/include/ops_stl.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

class DorofeevIRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    // 1. Достаем НАШ кортеж параметров из обертки фреймворка
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    
    // 2. И уже из нашего кортежа берем размер массива (нулевой элемент)
    int size = std::get<0>(params);

    // Настраиваем генератор случайных чисел
    std::mt19937 gen(42); // Фиксированный сид для стабильности тестов
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);

    input_data_.resize(size);
    for (int i = 0; i < size; ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Создаем эталонный вектор и сортируем его стандартным методом
    std::vector<double> expected = input_data_;
    std::sort(expected.begin(), expected.end());
    
    // Сравниваем результат твоей таски с эталоном
    return output_data == expected;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(DorofeevIRunFuncTestsThreads, TestSort) {
  ExecuteTest(GetParam());
}

// Задаем параметры: размер массива и строковое описание для логов
const std::array<TestType, 4> kTestParam = {
    std::make_tuple(10, "Small_Array"), 
    std::make_tuple(128, "Power_Of_Two"), 
    std::make_tuple(137, "Odd_Size"),
    std::make_tuple(1000, "Large_Array")
};

const auto kTaskName = PPC_SETTINGS_dorofeev_i_bitwise_sort_double_eo_batcher_merge;

// Собираем все реализации (ALL, OMP, SEQ, STL, TBB) в один тестовый набор
const auto kTestTasksList =
    std::tuple_cat(/*ppc::util::AddFuncTask<DorofeevIBitwiseSortDoubleEOBatcherMergeALL, InType>(kTestParam, kTaskName),
                   ppc::util::AddFuncTask<DorofeevIBitwiseSortDoubleEOBatcherMergeOMP, InType>(kTestParam, kTaskName),*/
                   ppc::util::AddFuncTask<DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ, InType>(kTestParam, kTaskName)
                   /*ppc::util::AddFuncTask<DorofeevIBitwiseSortDoubleEOBatcherMergeSTL, InType>(kTestParam, kTaskName),
                   ppc::util::AddFuncTask<DorofeevIBitwiseSortDoubleEOBatcherMergeTBB, InType>(kTestParam, kTaskName)*/);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = DorofeevIRunFuncTestsThreads::PrintFuncTestName<DorofeevIRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SortTests, DorofeevIRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge