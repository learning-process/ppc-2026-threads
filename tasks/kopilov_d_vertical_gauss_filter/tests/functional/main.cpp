#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "kopilov_d_vertical_gauss_filter/omp/include/ops_omp.hpp"
#include "kopilov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "kopilov_d_vertical_gauss_filter/stl/include/ops_stl.hpp"
#include "kopilov_d_vertical_gauss_filter/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kopilov_d_vertical_gauss_filter {

class VerticalGaussFilterTaskTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    const auto &matrix = std::get<0>(param);
    std::string res = std::to_string(matrix.width) + "x" + std::to_string(matrix.height);
    if (matrix.data.empty()) {
      res.append("_empty");
    } else {
      res.append("_v").append(std::to_string(matrix.data[0]));
    }
    return res;
  }

 protected:
  void SetUp() override {
    const auto &test_params = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_ = std::get<0>(test_params);
    reference_ = std::get<1>(test_params);
  }

  bool CheckTestOutputData(OutType &actual) final {
    return actual.width == reference_.width && actual.height == reference_.height && actual.data == reference_.data;
  }

  InType GetTestInputData() final {
    return input_;
  }

 private:
  InType input_;
  OutType reference_;
};

TEST_P(VerticalGaussFilterTaskTest, ValidationAndRun) {
  ExecuteTest(GetParam());
}

namespace {

Matrix CreateMatrix(int w, int h, const std::vector<uint8_t> &raw) {
  return Matrix{.width = w, .height = h, .data = raw};
}

const std::array<TestType, 5> kFunctionalScenarios = {
    std::make_tuple(CreateMatrix(1, 1, {100}), CreateMatrix(1, 1, {100})),
    std::make_tuple(CreateMatrix(2, 2, {1, 2, 3, 4}), CreateMatrix(2, 2, {1, 2, 2, 3})),
    std::make_tuple(CreateMatrix(3, 3, std::vector<uint8_t>(9, 16)), CreateMatrix(3, 3, std::vector<uint8_t>(9, 16))),
    std::make_tuple(CreateMatrix(3, 3, std::vector<uint8_t>(9, 42)), CreateMatrix(3, 3, std::vector<uint8_t>(9, 42))),
    std::make_tuple(CreateMatrix(4, 4, std::vector<uint8_t>(16, 100)),
                    CreateMatrix(4, 4, std::vector<uint8_t>(16, 100)))};

const auto kTaskProviders = std::tuple_cat(ppc::util::AddFuncTask<KopilovDVerticalGaussFilterSEQ, InType>(
                                               kFunctionalScenarios, PPC_SETTINGS_kopilov_d_vertical_gauss_filter),
                                           ppc::util::AddFuncTask<KopilovDVerticalGaussFilterOMP, InType>(
                                               kFunctionalScenarios, PPC_SETTINGS_kopilov_d_vertical_gauss_filter),
                                           ppc::util::AddFuncTask<KopilovDVerticalGaussFilterSTL, InType>(
                                               kFunctionalScenarios, PPC_SETTINGS_kopilov_d_vertical_gauss_filter),
                                           ppc::util::AddFuncTask<KopilovDVerticalGaussFilterTBB, InType>(
                                               kFunctionalScenarios, PPC_SETTINGS_kopilov_d_vertical_gauss_filter));

INSTANTIATE_TEST_SUITE_P(GaussFilters, VerticalGaussFilterTaskTest, ppc::util::ExpandToValues(kTaskProviders),
                         VerticalGaussFilterTaskTest::PrintFuncTestName<VerticalGaussFilterTaskTest>);

}  // namespace

namespace {

struct ValidationCase {
  std::function<std::shared_ptr<BaseTask>(const Matrix &)> generator;
  Matrix test_input;
  std::string tag;
};

void PrintTo(const ValidationCase &c, std::ostream *os) {
  *os << c.tag;
}

class VerticalGaussFilterValidationTest : public testing::TestWithParam<ValidationCase> {};

TEST_P(VerticalGaussFilterValidationTest, DetectsInvalidInput) {
  const auto &p = GetParam();
  EXPECT_FALSE(p.generator(p.test_input)->Validation());
}

std::vector<ValidationCase> GetValidationData() {
  std::vector<ValidationCase> cases;

  std::vector<std::pair<std::string, Matrix>> bad_inputs = {
      {"ZeroSize", {.width = 0, .height = 0, .data = {}}},
      {"ZeroWidth", {.width = 0, .height = 10, .data = std::vector<uint8_t>(10, 0)}},
      {"ZeroHeight", {.width = 10, .height = 0, .data = std::vector<uint8_t>(10, 0)}},
      {"SizeMismatch", {.width = 3, .height = 3, .data = {255, 128}}},
      {"NegativeW", {.width = -1, .height = 1, .data = {0}}},
      {"NegativeH", {.width = 1, .height = -1, .data = {0}}}};

  auto register_impl = [&](const std::string &name, auto factory) {
    for (const auto &item : bad_inputs) {
      std::string full_tag = name;
      full_tag.append("_").append(item.first);
      cases.push_back({factory, item.second, full_tag});
    }
  };

  register_impl("SEQ", [](const Matrix &m) { return std::make_shared<KopilovDVerticalGaussFilterSEQ>(m); });
  register_impl("OMP", [](const Matrix &m) { return std::make_shared<KopilovDVerticalGaussFilterOMP>(m); });
  register_impl("STL", [](const Matrix &m) { return std::make_shared<KopilovDVerticalGaussFilterSTL>(m); });
  register_impl("TBB", [](const Matrix &m) { return std::make_shared<KopilovDVerticalGaussFilterTBB>(m); });

  return cases;
}

INSTANTIATE_TEST_SUITE_P(NegativeTests, VerticalGaussFilterValidationTest, testing::ValuesIn(GetValidationData()),
                         [](const testing::TestParamInfo<ValidationCase> &info) { return info.param.tag; });

}  // namespace

}  // namespace kopilov_d_vertical_gauss_filter
