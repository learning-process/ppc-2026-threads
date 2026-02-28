#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "kazennova_a_fox_algorithm/common/include/common.hpp"
#include "kazennova_a_fox_algorithm/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kazennova_a_fox_algorithm {

class KazennovaAFuncTestSeq : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-10.0, 10.0);

    // Матрица A (size × size)
    input_data_.A.rows = size;
    input_data_.A.cols = size;
    input_data_.A.data.resize(static_cast<size_t>(size * size));
    for (int i = 0; i < size * size; ++i) {
      input_data_.A.data[i] = dis(gen);
    }

    // Матрица B (size × size)
    input_data_.B.rows = size;
    input_data_.B.cols = size;
    input_data_.B.data.resize(static_cast<size_t>(size * size));
    for (int i = 0; i < size * size; ++i) {
      input_data_.B.data[i] = dis(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Проверяем размер результирующей матрицы
    if (output_data.rows != input_data_.A.rows || output_data.cols != input_data_.B.cols) {
      return false;
    }
    if (output_data.data.size() != static_cast<size_t>(output_data.rows * output_data.cols)) {
      return false;
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

TEST_P(KazennovaAFuncTestSeq, MatrixMultiplicationTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {std::make_tuple(2, "2x2"), std::make_tuple(3, "3x3"),
                                            std::make_tuple(5, "5x5"), std::make_tuple(10, "10x10")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KazennovaATestTaskSEQ, InType>(kTestParam, PPC_SETTINGS_kazennova_a_fox_algorithm));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KazennovaAFuncTestSeq::PrintFuncTestName<KazennovaAFuncTestSeq>;

INSTANTIATE_TEST_SUITE_P(FoxSeqTests, KazennovaAFuncTestSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kazennova_a_fox_algorithm
