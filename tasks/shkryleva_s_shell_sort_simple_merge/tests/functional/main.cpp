#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"
#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/omp/include/ops_omp.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "shkryleva_s_shell_sort_simple_merge/stl/include/ops_stl.hpp"
#include "shkryleva_s_shell_sort_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergeFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 protected:
  void SetUp() override {
    TestType param = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    input_data_ = std::get<0>(param);
    expected_data_ = std::get<1>(param);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_data_.size()) {
      return false;
    }

    for (std::size_t i = 0; i < output_data.size(); i++) {
      if (output_data[i] != expected_data_[i]) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 public:
  static std::string PrintTestParam(const TestType &param) {
    const auto &input = std::get<0>(param);
    return "Size_" + std::to_string(input.size());
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

}  // namespace shkryleva_s_shell_sort_simple_merge

namespace {

using namespace shkryleva_s_shell_sort_simple_merge;

TEST_P(ShkrylevaSShellMergeFuncTests, shellMergeTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    TestType{InType{1, 2, 3, 4, 5, 6}, OutType{1, 2, 3, 4, 5, 6}},
    TestType{InType{9, 8, 7, 6, 5, 4, 3, 2, 1}, OutType{1, 2, 3, 4, 5, 6, 7, 8, 9}},
    TestType{InType{1}, OutType{1}},
    TestType{InType{2, 2, 2, 2, 2, 2, 2, 2}, OutType{2, 2, 2, 2, 2, 2, 2, 2}},
    TestType{InType{2, 2, 44, 2, 3, 5, 1}, OutType{1, 2, 2, 2, 3, 5, 44}},
    TestType{InType{2, 1}, OutType{1, 2}},
    TestType{InType{1, -2, 3, -5}, OutType{-5, -2, 1, 3}},
    TestType{InType{1, 22, 13, 51, 2, 1, 2, 2, 34, 41}, OutType{1, 1, 2, 2, 2, 13, 22, 34, 41, 51}}};

// Создаём списки задач для каждой версии
const auto kTestTasksListSeq = ppc::util::AddFuncTask<ShkrylevaSShellMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kTestTasksListOmp = ppc::util::AddFuncTask<ShkrylevaSShellMergeOMP, InType>(
    kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kTestTasksListStl = ppc::util::AddFuncTask<ShkrylevaSShellMergeSTL, InType>(
    kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kTestTasksListTbb = ppc::util::AddFuncTask<ShkrylevaSShellMergeTBB, InType>(
    kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kTestTasksListAll = ppc::util::AddFuncTask<ShkrylevaSShellMergeALL, InType>(
    kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);

// Объединяем в один вектор
std::vector<ppc::task::Task<InType, OutType>> kAllTestTasks;
kAllTestTasks.reserve(kTestTasksListSeq.size() + kTestTasksListOmp.size() + kTestTasksListStl.size() +
                      kTestTasksListTbb.size() + kTestTasksListAll.size());
kAllTestTasks.insert(kAllTestTasks.end(), kTestTasksListSeq.begin(), kTestTasksListSeq.end());
kAllTestTasks.insert(kAllTestTasks.end(), kTestTasksListOmp.begin(), kTestTasksListOmp.end());
kAllTestTasks.insert(kAllTestTasks.end(), kTestTasksListStl.begin(), kTestTasksListStl.end());
kAllTestTasks.insert(kAllTestTasks.end(), kTestTasksListTbb.begin(), kTestTasksListTbb.end());
kAllTestTasks.insert(kAllTestTasks.end(), kTestTasksListAll.begin(), kTestTasksListAll.end());

const auto kGtestValues = ppc::util::ExpandToValues(kAllTestTasks);
const auto kTestName = ShkrylevaSShellMergeFuncTests::PrintFuncTestName<ShkrylevaSShellMergeFuncTests>;

INSTANTIATE_TEST_SUITE_P(shellMergeTests, ShkrylevaSShellMergeFuncTests, kGtestValues, kTestName);

}  // namespace
