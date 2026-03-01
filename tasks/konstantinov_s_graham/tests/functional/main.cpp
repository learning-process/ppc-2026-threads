#include <gtest/gtest.h>
#include <stb/stb_image.h>

// #include <algorithm>
#include <array>
#include <cstddef>
// #include <cstdint>
// #include <numeric>
// #include <stdexcept>
#include <string>
#include <tuple>
// #include <utility>
#include <vector>

// #include "konstantinov_s_graham/all/include/ops_all.hpp"
#include "konstantinov_s_graham/common/include/common.hpp"
// #include "konstantinov_s_graham/omp/include/ops_omp.hpp"
#include "konstantinov_s_graham/seq/include/ops_seq.hpp"
// #include "konstantinov_s_graham/stl/include/ops_stl.hpp"
// #include "konstantinov_s_graham/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace konstantinov_s_graham {

class KonstantinovSRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param).first.size()) + "_" + std::to_string(std::get<1>(test_param).size());
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    test_input_ = std::get<0>(params);
    test_expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // std::cout<<std::endl;
    // if (output_data.size() != test_expected_output_.size()) {
    //   // std::cout<<"AAAAA "<<output_data.size()<<"out< >ans "<<test_expected_output_.size()<<std::endl;
    //   // return false;
    // }
    // for (int i = 0; (i < std::max(output_data.size(), test_expected_output_.size())) && (output_data.size() > 0);
    // i++) {
    //   // if (i < output_data.size()) {
    //   //  std::cout<<output_data[i].first << " " << output_data[i].second;
    //   //  std::cout<< " -- ";
    //   // if (test_expected_output_.size()) }
    //   // std::cout << test_expected_output_[i].first<< " "<<test_expected_output_[i].second;
    //   // std::cout<<std::endl;
    // }

    return test_expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return test_input_;
  }

 private:
  InType test_input_;
  OutType test_expected_output_;
};

namespace {

TEST_P(KonstantinovSRunFuncTestsThreads, GrahamTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {

    std::make_tuple(InType{{1.0}, {1.0}}, OutType{{1.0, 1.0}}),

    std::make_tuple(InType{{}, {}}, OutType{}),

    std::make_tuple(InType{{1.0, 0.0}, {1.0, 0.0}}, OutType{{0.0, 0.0}, {1.0, 1.0}}),

    std::make_tuple(InType{{0.0, 1.0, 2.0, 1.0}, {0.0, 1.0, 2.0, 1.0}}, OutType{{0.0, 0.0}, {2.0, 2.0}}),

    std::make_tuple(InType{{0.0, 0.0, 1.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 1.0}},
                    OutType{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}),

    std::make_tuple(InType{{0.0, 2.0, 1.0, -1.0, -2.0, 0.0, 1.0}, {0.0, 1.0, -1.0, -1.0, 2.0, 3.0, 1.0}},
                    OutType{{-1.0, -1.0}, {1.0, -1.0}, {2.0, 1.0}, {0.0, 3.0}, {-2.0, 2.0}}),

    std::make_tuple(InType{{5.5, -2.1, 3.3, -4.4, 0.0, -1.5}, {-3.2, 4.8, 3.3, -4.4, 2.2, 0.0}},
                    OutType{{-4.4, -4.4}, {5.5, -3.2}, {3.3, 3.3}, {-2.1, 4.8}}),

    std::make_tuple(InType{{0.1, -0.2, 0.4, 1.1, 0.0}, {0.5, 0.3, -0.6, 1.1, -1.0}},
                    OutType{{0.0, -1.0}, {0.4, -0.6}, {1.1, 1.1}, {0.1, 0.5}, {-0.2, 0.3}}),

    std::make_tuple(InType{{9e9, -9e9, 9e9, -9e9, 0.0}, {9e9, 9e9, -9e9, -9e9, 0.0}},
                    OutType{{-9e9, -9e9}, {9e9, -9e9}, {9e9, 9e9}, {-9e9, 9e9}}),

    std::make_tuple(
        InType{{0.0, 1.0, 1.0, 0.0, 0.2, 0.3, 0.7, 0.6, 0.4}, {0.0, 0.0, 1.0, 1.0, 0.2, 0.7, 0.3, 0.6, 0.5}},
        OutType{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}})};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KonstantinovAGrahamSEQ, InType>(kTestParam, PPC_SETTINGS_konstantinov_s_graham));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KonstantinovSRunFuncTestsThreads::PrintFuncTestName<KonstantinovSRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(GrahamTests, KonstantinovSRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace konstantinov_s_graham
