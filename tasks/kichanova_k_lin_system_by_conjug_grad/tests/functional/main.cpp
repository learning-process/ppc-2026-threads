#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"
#include "kichanova_k_lin_system_by_conjug_grad/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

namespace {

LinSystemData CreateIdentitySystem(int n, size_t n_squared) {
  LinSystemData data;
  data.n = n;
  data.epsilon = 1e-10;
  data.A.assign(n_squared, 0.0);
  const auto stride = static_cast<size_t>(n);
  for (int i = 0; i < n; ++i) {
    data.A[(static_cast<size_t>(i) * stride) + i] = 1.0;
  }
  data.b.assign(static_cast<size_t>(n), 1.0);
  return data;
}

LinSystemData CreateDiagonalSystem(int n, size_t n_squared) {
  LinSystemData data;
  data.n = n;
  data.epsilon = 1e-10;
  data.A.assign(n_squared, 0.0);
  const auto stride = static_cast<size_t>(n);
  for (int i = 0; i < n; ++i) {
    data.A[(static_cast<size_t>(i) * stride) + i] = static_cast<double>(i + 1);
  }
  data.b.assign(static_cast<size_t>(n), 1.0);
  return data;
}

LinSystemData CreateTridiagonalSystem(int n, size_t n_squared) {
  LinSystemData data;
  data.n = n;
  data.epsilon = 1e-10;
  data.A.assign(n_squared, 0.0);
  const auto stride = static_cast<size_t>(n);
  for (int i = 0; i < n; ++i) {
    data.A[(static_cast<size_t>(i) * stride) + i] = 2.0;
    if (i > 0) {
      data.A[(static_cast<size_t>(i) * stride) + (i - 1)] = -1.0;
    }
    if (i < n - 1) {
      data.A[(static_cast<size_t>(i) * stride) + (i + 1)] = -1.0;
    }
  }
  data.b.resize(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    data.b[i] = static_cast<double>(i + 1);
  }
  return data;
}

LinSystemData CreateRandomSPDSystem(int n, size_t n_squared) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1.0, 1.0);

  LinSystemData data;
  data.n = n;
  data.epsilon = 1e-10;

  std::vector<double> m(n_squared);
  for (int i = 0; i < n * n; ++i) {
    m[i] = dis(gen);
  }

  data.A.assign(n_squared, 0.0);
  const auto stride = static_cast<size_t>(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += m[(static_cast<size_t>(i) * stride) + k] * m[(static_cast<size_t>(j) * stride) + k];
      }
      data.A[(static_cast<size_t>(i) * stride) + j] = sum + n;
    }
  }

  data.b.resize(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    data.b[i] = dis(gen);
  }

  return data;
}

LinSystemData CreateTestSystem(int n, const std::string &type) {
  const size_t n_squared = static_cast<size_t>(n) * n;

  if (type == "identity") {
    return CreateIdentitySystem(n, n_squared);
  }
  if (type == "diagonal") {
    return CreateDiagonalSystem(n, n_squared);
  }
  if (type == "tridiagonal") {
    return CreateTridiagonalSystem(n, n_squared);
  }
  return CreateRandomSPDSystem(n, n_squared);
}

}  // namespace

class KichanovaKRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<LinSystemData, std::vector<double>, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "n" + std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int n = std::get<0>(params);
    std::string type = std::get<1>(params);

    input_data_ = CreateTestSystem(n, type);
  }

  bool CheckTestOutputData(std::vector<double> &output_data) final {
    if (output_data.size() != static_cast<size_t>(input_data_.n)) {
      return false;
    }

    double residual_norm = 0.0;
    const auto stride = static_cast<size_t>(input_data_.n);
    for (int i = 0; i < input_data_.n; ++i) {
      double sum = 0.0;
      for (int j = 0; j < input_data_.n; ++j) {
        const size_t pos = (static_cast<size_t>(i) * stride) + j;
        sum += input_data_.A[pos] * output_data[j];
      }
      double diff = sum - input_data_.b[i];
      residual_norm += diff * diff;
    }
    residual_norm = std::sqrt(residual_norm);

    return residual_norm < input_data_.epsilon * std::sqrt(static_cast<double>(input_data_.n));
  }

  LinSystemData GetTestInputData() final {
    return input_data_;
  }

 private:
  LinSystemData input_data_;
};

namespace {

TEST_P(KichanovaKRunFuncTestsThreads, SolveLinearSystem) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {
    std::make_tuple(2, "identity"),    std::make_tuple(3, "identity"),    std::make_tuple(5, "identity"),
    std::make_tuple(2, "diagonal"),    std::make_tuple(4, "diagonal"),    std::make_tuple(6, "diagonal"),
    std::make_tuple(3, "tridiagonal"), std::make_tuple(5, "tridiagonal"), std::make_tuple(7, "tridiagonal"),
    std::make_tuple(4, "random_spd"),  std::make_tuple(8, "random_spd"),  std::make_tuple(10, "random_spd")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<KichanovaKLinSystemByConjugGradSEQ, LinSystemData>(
    kTestParam, PPC_SETTINGS_kichanova_k_lin_system_by_conjug_grad));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = KichanovaKRunFuncTestsThreads::PrintFuncTestName<KichanovaKRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(LinearSystemTests, KichanovaKRunFuncTestsThreads, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace kichanova_k_lin_system_by_conjug_grad
