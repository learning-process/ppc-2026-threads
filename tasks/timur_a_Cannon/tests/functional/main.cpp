#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "timur_a_Cannon/common/include/common.hpp"
#include "timur_a_Cannon/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace timur_a_Cannon {

namespace {

Matrix NaiveMultiply(const Matrix &A, const Matrix &B) {
  Matrix C(A.n, 0.0);
  for (std::size_t i = 0; i < A.n; ++i) {
    for (std::size_t k = 0; k < A.n; ++k) {
      double a_val = A(i, k);
      for (std::size_t j = 0; j < A.n; ++j) {
        C(i, j) += a_val * B(k, j);
      }
    }
  }
  return C;
}

TaskData MakeTestMatrices(std::size_t n, const std::string &pattern) {
  TaskData data;
  data.A = Matrix(n, 0.0);
  data.B = Matrix(n, 0.0);

  if (pattern == "identity") {
    for (std::size_t i = 0; i < n; ++i) {
      data.A(i, i) = 1.0;
      data.B(i, i) = 1.0;
    }
  } else if (pattern == "ones") {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        data.A(i, j) = 1.0;
        data.B(i, j) = 1.0;
      }
    }
  } else {
    double val = 1.0;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        data.A(i, j) = val;
        data.B(i, j) = val + 1.0;
        val += 1.0;
      }
    }
  }

  return data;
}

}  // namespace

class TimurACannonFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    auto n = std::get<0>(params);
    auto pattern = std::get<1>(params);

    input_data_ = MakeTestMatrices(n, pattern);
    reference_ = NaiveMultiply(input_data_.A, input_data_.B);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == reference_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType reference_{};
};

namespace {

TEST_P(TimurACannonFuncTests, MatrixMultiplyCorrectness) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(2, "identity"), std::make_tuple(4, "ones"),
                                            std::make_tuple(6, "sequence")};

const auto kTestTasksList = ppc::util::AddFuncTask<TimurACannonSEQ, InType>(kTestParam, PPC_SETTINGS_timur_a_Cannon);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = TimurACannonFuncTests::PrintTestParam;

INSTANTIATE_TEST_SUITE_P(CannonSeqTests, TimurACannonFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace timur_a_Cannon
