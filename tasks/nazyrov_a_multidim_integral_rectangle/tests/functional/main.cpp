#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "nazyrov_a_multidim_integral_rectangle/common/include/common.hpp"
#include "nazyrov_a_multidim_integral_rectangle/omp/include/ops_omp.hpp"
#include "nazyrov_a_multidim_integral_rectangle/seq/include/ops_seq.hpp"
#include "nazyrov_a_multidim_integral_rectangle/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nazyrov_a_multidim_integral_rectangle {
namespace {

constexpr double kTol = 1e-4;

struct TestDef {
  Func func;
  Bounds bounds;
  double expected{};
};

TestDef MakeTestDef(int func_id, int dim) {
  switch (func_id) {
    case 0: {
      // f(x) = 1, bounds [0,1]^d, expected = 1
      return {.func = [](const std::vector<double> &) {
        return 1.0;
      }, .bounds = Bounds(static_cast<std::size_t>(dim), {0.0, 1.0}), .expected = 1.0};
    }
    case 1: {
      // f(x) = sum(x_i), bounds [0,1]^d, expected = d/2
      const double expected = 0.5 * static_cast<double>(dim);
      return {.func = [](const std::vector<double> &x) {
        double s = 0.0;
        for (double xi : x) {
          s += xi;
        }
        return s;
      }, .bounds = Bounds(static_cast<std::size_t>(dim), {0.0, 1.0}), .expected = expected};
    }
    case 2: {
      // f(x) = sum(x_i^2), bounds [0,1]^d, expected = d/3
      const double expected = static_cast<double>(dim) / 3.0;
      return {.func = [](const std::vector<double> &x) {
        double s = 0.0;
        for (double xi : x) {
          s += xi * xi;
        }
        return s;
      }, .bounds = Bounds(static_cast<std::size_t>(dim), {0.0, 1.0}), .expected = expected};
    }
    case 3: {
      // f(x) = prod(x_i), bounds [1,2]^d, expected = (3/2)^d
      const double expected = std::pow(1.5, static_cast<double>(dim));
      return {.func = [](const std::vector<double> &x) {
        double p = 1.0;
        for (double xi : x) {
          p *= xi;
        }
        return p;
      }, .bounds = Bounds(static_cast<std::size_t>(dim), {1.0, 2.0}), .expected = expected};
    }
    default:
      return {.func = [](const std::vector<double> &) {
        return 0.0;
      }, .bounds = Bounds(static_cast<std::size_t>(dim), {0.0, 1.0}), .expected = 0.0};
  }
}

}  // namespace

class NazyrovAMultidimIntegralRectangleFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &p) {
    return "f" + std::to_string(std::get<0>(p)) + "_d" + std::to_string(std::get<1>(p)) + "_n" +
           std::to_string(std::get<2>(p));
  }

 protected:
  void SetUp() override {
    const auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int func_id = std::get<0>(params);
    const int dim = std::get<1>(params);
    const int steps = std::get<2>(params);

    auto def = MakeTestDef(func_id, dim);
    expected_ = def.expected;
    input_data_ = std::make_tuple(def.func, def.bounds, steps);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    if (ppc::util::IsUnderMpirun()) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    if (rank != 0) {
      return true;
    }
    return std::abs(output_data - expected_) < kTol;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  double expected_{};
};

namespace {

TEST_P(NazyrovAMultidimIntegralRectangleFuncTests, ComputeMultidimIntegral) {
  ExecuteTest(GetParam());
}

// (func_id, dim, steps)
const std::array<TestType, 10> kTestParams = {
    std::make_tuple(0, 1, 100),  // constant, 1D
    std::make_tuple(0, 3, 50),   // constant, 3D
    std::make_tuple(1, 1, 100),  // linear sum, 1D
    std::make_tuple(1, 2, 100),  // linear sum, 2D
    std::make_tuple(1, 3, 50),   // linear sum, 3D
    std::make_tuple(2, 1, 200),  // quadratic, 1D
    std::make_tuple(2, 2, 200),  // quadratic, 2D
    std::make_tuple(2, 3, 100),  // quadratic, 3D
    std::make_tuple(3, 2, 100),  // product, 2D over [1,2]^2
    std::make_tuple(3, 3, 50),   // product, 3D over [1,2]^3
};

const auto kTaskList = std::tuple_cat(ppc::util::AddFuncTask<NazyrovAMultidimIntegralRectangleSeq, InType>(
                                          kTestParams, PPC_SETTINGS_nazyrov_a_multidim_integral_rectangle),
                                      ppc::util::AddFuncTask<NazyrovAMultidimIntegralRectangleOmp, InType>(
                                          kTestParams, PPC_SETTINGS_nazyrov_a_multidim_integral_rectangle),
                                      ppc::util::AddFuncTask<NazyrovAMultidimIntegralRectangleTbb, InType>(
                                          kTestParams, PPC_SETTINGS_nazyrov_a_multidim_integral_rectangle));

const auto kGtestValues = ppc::util::ExpandToValues(kTaskList);
const auto kPerfTestName =
    NazyrovAMultidimIntegralRectangleFuncTests::PrintFuncTestName<NazyrovAMultidimIntegralRectangleFuncTests>;

INSTANTIATE_TEST_SUITE_P(FuncTests, NazyrovAMultidimIntegralRectangleFuncTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace nazyrov_a_multidim_integral_rectangle
