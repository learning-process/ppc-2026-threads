#include <gtest/gtest.h>

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
#include <random>
#include <cmath>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"
#include "kichanova_k_lin_system_by_conjug_grad/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

inline LinSystemData CreateTestSystem(int n, const std::string& type) {
    LinSystemData data;
    data.n = n;
    data.epsilon = 1e-10;
    
    if (type == "identity") {
        data.A.assign(n * n, 0.0);
        for (int i = 0; i < n; ++i) {
            data.A[i * n + i] = 1.0;
        }
        data.b.assign(n, 1.0);
    }
    else if (type == "diagonal") {
        data.A.assign(n * n, 0.0);
        for (int i = 0; i < n; ++i) {
            data.A[i * n + i] = static_cast<double>(i + 1);
        }
        data.b.assign(n, 1.0);
    }
    else if (type == "tridiagonal") {
        data.A.assign(n * n, 0.0);
        for (int i = 0; i < n; ++i) {
            data.A[i * n + i] = 2.0;
            if (i > 0) data.A[i * n + (i - 1)] = -1.0;
            if (i < n - 1) data.A[i * n + (i + 1)] = -1.0;
        }
        data.b.resize(n);
        for (int i = 0; i < n; ++i) {
            data.b[i] = static_cast<double>(i + 1);
        }
    }
    else if (type == "random_spd") {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1.0, 1.0);
        
        std::vector<double> M(n * n);
        for (int i = 0; i < n * n; ++i) {
            M[i] = dis(gen);
        }
        
        data.A.assign(n * n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double sum = 0.0;
                for (int k = 0; k < n; ++k) {
                    sum += M[i * n + k] * M[j * n + k];
                }
                data.A[i * n + j] = sum + n;
            }
        }
        
        data.b.resize(n);
        for (int i = 0; i < n; ++i) {
            data.b[i] = dis(gen);
        }
    }
    
    return data;
}

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
    for (int i = 0; i < input_data_.n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < input_data_.n; ++j) {
            sum += input_data_.A[i * input_data_.n + j] * output_data[j];
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
    std::make_tuple(2, "identity"),
    std::make_tuple(3, "identity"),
    std::make_tuple(5, "identity"),
    std::make_tuple(2, "diagonal"),
    std::make_tuple(4, "diagonal"),
    std::make_tuple(6, "diagonal"),
    std::make_tuple(3, "tridiagonal"),
    std::make_tuple(5, "tridiagonal"),
    std::make_tuple(7, "tridiagonal"),
    std::make_tuple(4, "random_spd"),
    std::make_tuple(8, "random_spd"),
    std::make_tuple(10, "random_spd")
};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<KichanovaKLinSystemByConjugGradSEQ, LinSystemData>(kTestParam, PPC_SETTINGS_kichanova_k_lin_system_by_conjug_grad));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = KichanovaKRunFuncTestsThreads::PrintFuncTestName<KichanovaKRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(LinearSystemTests, KichanovaKRunFuncTestsThreads, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace kichanova_k_lin_system_by_conjug_grad