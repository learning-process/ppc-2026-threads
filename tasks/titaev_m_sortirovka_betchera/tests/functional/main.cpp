#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"
#include "titaev_m_sortirovka_betchera/common/include/common.hpp"
#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"
#include "titaev_m_sortirovka_betchera/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace titaev_m_sortirovka_betchera {

using ParamType =
    std::tuple<std::function<std::shared_ptr<ppc::task::Task<InType, OutType>>(const InType &)>, std::string, TestType>;

class TitaevBatcherRadixFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const testing::TestParamInfo<ParamType> &info) {
    return std::get<1>(info.param);
  }

 protected:
  InType input;
  void SetUp() override {
    ParamType full_param = GetParam();
    TestType param = std::get<2>(full_param);
    const int size = std::get<0>(param);
    std::mt19937_64 gen(size + 42);
    std::uniform_real_distribution<double> dist(-5000.0, 5000.0);
    input.resize(size);
    for (int i = 0; i < size; i++) {
      input[i] = dist(gen);
    }
  }
  bool CheckTestOutputData(OutType &output) final {
    if (output.size() != input.size()) {
      return false;
    }
    for (size_t i = 1; i < output.size(); i++) {
      if (output[i] < output[i - 1]) {
        return false;
      }
    }
    return true;
  }
  InType GetTestInputData() final {
    return input;
  }
};

namespace {
std::shared_ptr<ppc::task::Task<InType, OutType>> MakeSeqTask(const InType &in) {
  return std::make_shared<TitaevSortirovkaBetcheraSEQ>(in);
}
std::shared_ptr<ppc::task::Task<InType, OutType>> MakeOMPTask(const InType &in) {
  return std::make_shared<TitaevSortirovkaBetcheraOMP>(in);
}

const ParamType kParam1{MakeSeqTask, "seq_100", TestType{100, "100"}};
const ParamType kParam2{MakeSeqTask, "seq_500", TestType{500, "500"}};
const ParamType kParam3{MakeOMPTask, "omp_100", TestType{100, "100"}};
const ParamType kParam4{MakeOMPTask, "omp_500", TestType{500, "500"}};

INSTANTIATE_TEST_SUITE_P(FunctionalTests, TitaevBatcherRadixFuncTests,
                         ::testing::Values(kParam1, kParam2, kParam3, kParam4),
                         TitaevBatcherRadixFuncTests::PrintTestParam);
}  // namespace
}  // namespace titaev_m_sortirovka_betchera
