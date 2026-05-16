#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/seq/include/ops_seq.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

class VasilievMShellSortBatcherMergeFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    std::string name = std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
    for (char &c : name) {
      if (std::isalnum(static_cast<unsigned char>(c)) == 0) {
        c = '_';
      }
    }
    return name;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string filename = std::get<1>(params);

    std::string path =
        ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_vasiliev_m_shell_sort_batcher_merge), filename);

    std::ifstream file(path);

    if (!file.is_open()) {
      throw std::runtime_error("Wrong path.");
    }

    int size = 0;
    file >> size;

    std::vector<int> vec(size);

    for (int i = 0; i < size; i++) {
      file >> vec[i];
    }

    input_data_ = vec;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.empty() && !input_data_.empty()) {
      return false;
    }

    if (output_data.size() != input_data_.size()) {
      return false;
    }

    auto it = std::ranges::is_sorted_until(output_data);
    if (it != output_data.end()) {
      return false;
    }

    auto input_sorted = input_data_;
    std::ranges::sort(input_sorted);
    return input_sorted == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(VasilievMShellSortBatcherMergeFuncTestsSEQ, BatcherSorting) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParam = {std::make_tuple(0, "test_vec_1.txt"), std::make_tuple(1, "test_vec_2.txt"),
                                            std::make_tuple(2, "test_vec_3.txt"), std::make_tuple(3, "test_vec_4.txt"),
                                            std::make_tuple(4, "test_vec_5.txt"), std::make_tuple(5, "test_vec_6.txt")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<VasilievMShellSortBatcherMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_vasiliev_m_shell_sort_batcher_merge));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    VasilievMShellSortBatcherMergeFuncTestsSEQ::PrintFuncTestName<VasilievMShellSortBatcherMergeFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(FuncBatcherTests, VasilievMShellSortBatcherMergeFuncTestsSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace vasiliev_m_shell_sort_batcher_merge
