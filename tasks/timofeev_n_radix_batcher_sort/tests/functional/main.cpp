#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"
#include "timofeev_n_radix_batcher_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

class TimofeevNRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<3>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    for (int i = 0; std::cmp_less(i, static_cast<int>(output_data.size() - 1)); i++) {
      if (output_data[i] > output_data[i + 1]) {
        return false;
      }
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

TEST_P(TimofeevNRunFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

std::vector<int> t1 = {4, 3, 5};
std::vector<int> t2 = {1, -2, 3, -4, 5};
std::vector<int> t3 = {1, -2, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 13, -14, 15, -16};
std::vector<int> mt;

const std::array<TestType, 3> kTestParam = {std::make_tuple(t1, mt, 1, "3"), std::make_tuple(t2, mt, 1, "5"),
                                            std::make_tuple(t3, mt, 1, "16")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<TimofeevNRadixBatcherSEQ, InType>(kTestParam, PPC_SETTINGS_timofeev_n_radix_batcher_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = TimofeevNRunFuncTestsThreads::PrintFuncTestName<TimofeevNRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, TimofeevNRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace timofeev_n_radix_batcher_sort_threads
