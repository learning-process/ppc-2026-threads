#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "makovskiy_i_graham_hull/common/include/common.hpp"
#include "makovskiy_i_graham_hull/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace makovskiy_i_graham_hull {

class MakovskiyIGrahamHullRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);

    if (test_id == 1) {
      input_data_ = {{0, 0}, {0, 2}, {2, 0}, {2, 2}, {1, 1}};
    } else if (test_id == 2) {
      input_data_ = {{0, 0}, {5, 0}, {2, 5}};
    } else if (test_id == 3) {
      input_data_ = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {0, 3}, {3, 0}};
    } else if (test_id == 4) {
      input_data_ = {};
    } else if (test_id == 5) {
      input_data_ = {{1, 1}, {2, 2}};
    } else if (test_id == 6) {
      input_data_ = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}};
    } else if (test_id == 7) {
      input_data_ = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}};
    } else if (test_id == 8) {
      input_data_ = {{-2, -2}, {2, -2}, {2, 2}, {-2, 2}, {0, 0}};
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int test_id = std::get<0>(std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));

    if (test_id == 1) {
      return output_data.size() == 4;
    }
    if (test_id == 2) {
      return output_data.size() == 3;
    }
    if (test_id == 3) {
      return output_data.size() == 4;
    }
    if (test_id == 4) {
      return output_data.empty();
    }
    if (test_id == 5) {
      return output_data.size() == 2;
    }
    if (test_id == 6) {
      return output_data.size() == 2;
    }
    if (test_id == 7) {
      return output_data.size() == 2;
    }
    if (test_id == 8) {
      return output_data.size() == 4;
    }

    return false;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(MakovskiyIGrahamHullRunFuncTestsThreads, GrahamHullTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {std::make_tuple(1, "square_with_internal_point"),
                                            std::make_tuple(2, "triangle"),
                                            std::make_tuple(3, "collinear_points_on_edges"),
                                            std::make_tuple(4, "empty_array"),
                                            std::make_tuple(5, "two_points"),
                                            std::make_tuple(6, "diagonal_line"),
                                            std::make_tuple(7, "vertical_line"),
                                            std::make_tuple(8, "negative_coordinates")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<ConvexHullGrahamSEQ, InType>(kTestParam, PPC_SETTINGS_makovskiy_i_graham_hull));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    MakovskiyIGrahamHullRunFuncTestsThreads::PrintFuncTestName<MakovskiyIGrahamHullRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(GrahamHullTestsSuite, MakovskiyIGrahamHullRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace makovskiy_i_graham_hull
