#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"
#include "kruglova_a_conjugate_gradient_sle/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace kruglova_a_conjugate_gradient_sle {

class KruglovaAFuncTestAConjGradSle : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param) + "_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);
    std::string type = std::get<1>(params);

    input_data_.size = size;
    input_data_.A.assign(static_cast<std::size_t>(size) * static_cast<std::size_t>(size), 0.0);
    input_data_.b.assign(size, 0.0);

    if (type == "Identity") {
      FillIdentity(size);
    } else if (type == "Diagonal") {
      FillDiagonal(size);
    } else {
      FillSpd(size);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int n = input_data_.size;
    double max_err = 0.0;
    for (int i = 0; i < n; ++i) {
      double ax = 0.0;
      for (int j = 0; j < n; ++j) {
        ax += input_data_.A[(static_cast<std::size_t>(i) * n) + j] * output_data[j];
      }
      max_err = std::max(max_err, std::abs(ax - input_data_.b[i]));
    }
    return max_err < 1e-4;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;

  void FillIdentity(int size) {
    for (int i = 0; i < size; ++i) {
      input_data_.A[(static_cast<std::size_t>(i) * size) + i] = 1.0;
      input_data_.b[i] = static_cast<double>(i + 1);
    }
  }

  void FillDiagonal(int size) {
    for (int i = 0; i < size; ++i) {
      input_data_.A[(static_cast<std::size_t>(i) * size) + i] = static_cast<double>(i + 1) * 10.0;
      input_data_.b[i] = input_data_.A[(static_cast<std::size_t>(i) * size) + i];
    }
  }

  void FillSpd(int size) {
    for (int i = 0; i < size; ++i) {
      double sum = 0.0;
      for (int j = 0; j < size; ++j) {
        if (i != j) {
          auto val = static_cast<double>(((i + j) % 5) + 1);
          input_data_.A[(static_cast<std::size_t>(i) * size) + j] = val;
          input_data_.A[(static_cast<std::size_t>(j) * size) + i] = val;
          sum += val;
        }
      }
      input_data_.A[(static_cast<std::size_t>(i) * size) + i] = sum + 10.0;
      input_data_.b[i] = static_cast<double>((i % 3) + 1);
    }
  }
};

TEST_P(KruglovaAFuncTestAConjGradSle, SolveSystem) {
  ExecuteTest(GetParam());
}

namespace {
const std::array<TestType, 6> kTestParam = {std::make_tuple(1, "Size1"),    std::make_tuple(3, "Identity"),
                                            std::make_tuple(5, "Diagonal"), std::make_tuple(10, "Size10"),
                                            std::make_tuple(50, "Size50"),  std::make_tuple(101, "OddSize")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<KruglovaAConjGradSleSEQ, InType>(
    kTestParam, PPC_SETTINGS_kruglova_a_conjugate_gradient_sle));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = KruglovaAFuncTestAConjGradSle::PrintFuncTestName<KruglovaAFuncTestAConjGradSle>;

INSTANTIATE_TEST_SUITE_P(SleSolverTests, KruglovaAFuncTestAConjGradSle, kGtestValues, kPerfTestName);
}  // namespace

}  // namespace kruglova_a_conjugate_gradient_sle
