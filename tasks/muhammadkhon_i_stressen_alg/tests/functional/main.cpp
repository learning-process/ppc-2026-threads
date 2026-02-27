#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "muhammadkhon_i_stressen_alg/common/include/common.hpp"
#include "muhammadkhon_i_stressen_alg/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace muhammadkhon_i_stressen_alg {

class MuhammadkhonIStressenAlgFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int n = std::get<0>(params);
    const int size = n * n;

    std::vector<double> a(static_cast<size_t>(size));
    std::vector<double> b(static_cast<size_t>(size));
    for (int i = 0; i < size; ++i) {
      a[static_cast<size_t>(i)] = static_cast<double>((i % 7) + 1);
      b[static_cast<size_t>(i)] = static_cast<double>(((i * 3 + 5) % 11) + 1);
    }

    input_data_ = MatrixInput{.a = a, .b = b, .n = n};

    expected_output_.assign(static_cast<size_t>(size), 0.0);
    for (int row = 0; row < n; ++row) {
      for (int k = 0; k < n; ++k) {
        for (int col = 0; col < n; ++col) {
          expected_output_[static_cast<size_t>((static_cast<ptrdiff_t>(row) * n) + col)] +=
              a[static_cast<size_t>((static_cast<ptrdiff_t>(row) * n) + k)] *
              b[static_cast<size_t>((static_cast<ptrdiff_t>(k) * n) + col)];
        }
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }
    constexpr double kEps = 1e-9;
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::fabs(output_data[i] - expected_output_[i]) > kEps) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType expected_output_;
};

namespace {

TEST_P(MuhammadkhonIStressenAlgFuncTests, StrassenMatmul) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(2, "2"), std::make_tuple(3, "3"),   std::make_tuple(5, "5"),
    std::make_tuple(7, "7"), std::make_tuple(16, "16"),
};

const auto kTestTasksList =
    ppc::util::AddFuncTask<MuhammadkhonIStressenAlgSEQ, InType>(kTestParam, PPC_SETTINGS_muhammadkhon_i_stressen_alg);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = MuhammadkhonIStressenAlgFuncTests::PrintFuncTestName<MuhammadkhonIStressenAlgFuncTests>;

INSTANTIATE_TEST_SUITE_P(SeqMatrixTests, MuhammadkhonIStressenAlgFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace muhammadkhon_i_stressen_alg
