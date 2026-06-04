#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "titaev_m_sortirovka_betchera/all/include/ops_all.hpp"
#include "titaev_m_sortirovka_betchera/common/include/common.hpp"
#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"
#include "titaev_m_sortirovka_betchera/seq/include/ops_seq.hpp"
#include "titaev_m_sortirovka_betchera/stl/include/ops_stl.hpp"
#include "titaev_m_sortirovka_betchera/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace titaev_m_sortirovka_betchera {

class TitaevSortirovkaBetcheraFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    size_ = std::get<0>(params);
    input_data_.resize(static_cast<std::size_t>(size_));
    for (int i = 0; i < size_; i++) {
      input_data_[i] = static_cast<double>(size_ - i) - 0.5;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != static_cast<std::size_t>(size_)) {
      return false;
    }
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  int size_ = 0;
  InType input_data_;
};

namespace {
TEST_P(TitaevSortirovkaBetcheraFuncTests, SortCorrectness) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    std::make_tuple(1, "N1"),   std::make_tuple(2, "N2"),   std::make_tuple(4, "N4"), std::make_tuple(8, "N8"),
    std::make_tuple(16, "N16"), std::make_tuple(32, "N32"), std::make_tuple(7, "N7"), std::make_tuple(100, "N100")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<TitaevSortirovkaBetcheraSEQ, InType>(kTestParam, PPC_SETTINGS_titaev_m_sortirovka_betchera),
    ppc::util::AddFuncTask<TitaevSortirovkaBetcheraOMP, InType>(kTestParam, PPC_SETTINGS_titaev_m_sortirovka_betchera),
    ppc::util::AddFuncTask<TitaevSortirovkaBetcheraTBB, InType>(kTestParam, PPC_SETTINGS_titaev_m_sortirovka_betchera),
    ppc::util::AddFuncTask<TitaevSortirovkaBetcheraSTL, InType>(kTestParam, PPC_SETTINGS_titaev_m_sortirovka_betchera),
    ppc::util::AddFuncTask<TitaevSortirovkaBetcheraALL, InType>(kTestParam, PPC_SETTINGS_titaev_m_sortirovka_betchera));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kFuncTestName = TitaevSortirovkaBetcheraFuncTests::PrintFuncTestName<TitaevSortirovkaBetcheraFuncTests>;
INSTANTIATE_TEST_SUITE_P(SortFuncTests, TitaevSortirovkaBetcheraFuncTests, kGtestValues, kFuncTestName);
}  // namespace
}  // namespace titaev_m_sortirovka_betchera
