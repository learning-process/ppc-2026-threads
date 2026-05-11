#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <tuple>

#include "tochilin_e_hoar_sort_sim_mer/common/include/common.hpp"
#include "tochilin_e_hoar_sort_sim_mer/omp/include/ops_omp.hpp"
#include "tochilin_e_hoar_sort_sim_mer/seq/include/ops_seq.hpp"
#include "tochilin_e_hoar_sort_sim_mer/stl/include/ops_stl.hpp"
#include "tochilin_e_hoar_sort_sim_mer/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace tochilin_e_hoar_sort_sim_mer {

namespace {

std::uint32_t MakeSeed(int n, const std::string &desc) {
  std::uint32_t seed = 2166136261U;
  for (unsigned char ch : desc) {
    seed ^= ch;
    seed *= 16777619U;
  }
  seed ^= static_cast<std::uint32_t>(n);
  seed *= 16777619U;
  return seed;
}

}  // namespace

class TochilinEHoarSortSimMerRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    const int n = std::get<0>(params);
    const std::string &desc = std::get<1>(params);

    input_data_.resize(static_cast<std::size_t>(n));

    if (desc == "OneElement") {
      if (n > 0) {
        input_data_[0] = 42;
      }
    } else if (desc == "AlreadySorted") {
      for (int i = 0; i < n; ++i) {
        input_data_[static_cast<std::size_t>(i)] = i;
      }
    } else if (desc == "ReverseSorted") {
      for (int i = 0; i < n; ++i) {
        input_data_[static_cast<std::size_t>(i)] = n - i;
      }
    } else {
      if (n > 0) {
        std::mt19937 gen(MakeSeed(n, desc));
        std::uniform_int_distribution<> dis(-1000, 1000);
        for (int i = 0; i < n; ++i) {
          input_data_[static_cast<std::size_t>(i)] = dis(gen);
        }
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    OutType reference = input_data_;
    std::ranges::sort(reference);
    return reference == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;  // reopen
};

namespace {

TEST_P(TochilinEHoarSortSimMerRunFuncTests, TestSorting) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {
    std::make_tuple(1, "OneElement"),      std::make_tuple(2, "TwoElements"),  std::make_tuple(8, "EightElements"),
    std::make_tuple(13, "RandomSize"),     std::make_tuple(100, "MediumSize"), std::make_tuple(128, "AlreadySorted"),
    std::make_tuple(127, "ReverseSorted"), std::make_tuple(512, "LargeSize"),  std::make_tuple(1000, "LongSize"),
    std::make_tuple(2000, "LongLongSize")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<TochilinEHoarSortSimMerSTL, InType>(kTestParam, PPC_SETTINGS_tochilin_e_hoar_sort_sim_mer),
    ppc::util::AddFuncTask<TochilinEHoarSortSimMerTBB, InType>(kTestParam, PPC_SETTINGS_tochilin_e_hoar_sort_sim_mer),
    ppc::util::AddFuncTask<TochilinEHoarSortSimMerOMP, InType>(kTestParam, PPC_SETTINGS_tochilin_e_hoar_sort_sim_mer),
    ppc::util::AddFuncTask<TochilinEHoarSortSimMerSEQ, InType>(kTestParam, PPC_SETTINGS_tochilin_e_hoar_sort_sim_mer));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = TochilinEHoarSortSimMerRunFuncTests::PrintFuncTestName<TochilinEHoarSortSimMerRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(HoarSortSimMerTests, TochilinEHoarSortSimMerRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tochilin_e_hoar_sort_sim_mer
