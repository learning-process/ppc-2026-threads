#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "peterson_r_graham_scan/all/include/ops_all.hpp"
#include "peterson_r_graham_scan/common/include/common.hpp"
#include "peterson_r_graham_scan/omp/include/ops_omp.hpp"
#include "peterson_r_graham_scan/seq/include/ops_seq.hpp"
#include "peterson_r_graham_scan/stl/include/ops_stl.hpp"
#include "peterson_r_graham_scan/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace peterson_r_graham_scan {

class PetersonGrahamScannerFuncTests : public ppc::util::BaseRunFuncTests<InputValue, OutputValue, TestParameters> {
 public:
  static std::string PrintTestParam(const TestParameters &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestParameters params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutputValue &output_data) final {
    return input_data_ == output_data;
  }

  InputValue GetTestInputData() final {
    return input_data_;
  }

 private:
  InputValue input_data_ = 0;
};

namespace {

using Point = Point2D;

TEST_P(PetersonGrahamScannerFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestParameters, 3> kTestCases = {std::make_tuple(3, "circle_3"), std::make_tuple(5, "circle_5"),
                                                  std::make_tuple(7, "circle_7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<PetersonGrahamScannerSeq, InputValue>(kTestCases, PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::AddFuncTask<PetersonGrahamScannerOmp, InputValue>(kTestCases, PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::AddFuncTask<PetersonGrahamScannerTbb, InputValue>(kTestCases, PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::AddFuncTask<PetersonGrahamScannerStl, InputValue>(kTestCases, PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::AddFuncTask<PetersonGrahamScannerAll, InputValue>(kTestCases, PPC_SETTINGS_peterson_r_graham_scan));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestNameGenerator = PetersonGrahamScannerFuncTests::PrintFuncTestName<PetersonGrahamScannerFuncTests>;

INSTANTIATE_TEST_SUITE_P(DefaultTests, PetersonGrahamScannerFuncTests, kGtestValues, kTestNameGenerator);

void ExecutePipelineSeq(const std::shared_ptr<PetersonGrahamScannerSeq> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

void ExecutePipelineOmp(const std::shared_ptr<PetersonGrahamScannerOmp> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

void ExecutePipelineTbb(const std::shared_ptr<PetersonGrahamScannerTbb> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

void ExecutePipelineStl(const std::shared_ptr<PetersonGrahamScannerStl> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

void ExecutePipelineAll(const std::shared_ptr<PetersonGrahamScannerAll> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(PetersonGrahamScannerSeq, EmptyInput) {
  auto task = std::make_shared<PetersonGrahamScannerSeq>(0);
  ExecutePipelineSeq(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(PetersonGrahamScannerSeq, SinglePoint) {
  std::vector<Point> pts = {Point(5.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerSeq>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineSeq(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerSeq, AllIdenticalPoints) {
  std::vector<Point> pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerSeq>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineSeq(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerOmp, EmptyInput) {
  auto task = std::make_shared<PetersonGrahamScannerOmp>(0);
  ExecutePipelineOmp(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(PetersonGrahamScannerOmp, SinglePoint) {
  std::vector<Point> pts = {Point(5.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerOmp>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineOmp(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerOmp, AllIdenticalPoints) {
  std::vector<Point> pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerOmp>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineOmp(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerTbb, EmptyInput) {
  auto task = std::make_shared<PetersonGrahamScannerTbb>(0);
  ExecutePipelineTbb(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(PetersonGrahamScannerTbb, SinglePoint) {
  std::vector<Point> pts = {Point(5.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerTbb>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineTbb(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerTbb, AllIdenticalPoints) {
  std::vector<Point> pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerTbb>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineTbb(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerStl, EmptyInput) {
  auto task = std::make_shared<PetersonGrahamScannerStl>(0);
  ExecutePipelineStl(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(PetersonGrahamScannerStl, SinglePoint) {
  std::vector<Point> pts = {Point(5.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerStl>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineStl(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerStl, AllIdenticalPoints) {
  std::vector<Point> pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerStl>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineStl(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerAll, EmptyInput) {
  auto task = std::make_shared<PetersonGrahamScannerAll>(0);
  ExecutePipelineAll(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(PetersonGrahamScannerAll, SinglePoint) {
  std::vector<Point> pts = {Point(5.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerAll>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineAll(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerAll, AllIdenticalPoints) {
  std::vector<Point> pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerAll>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipelineAll(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

}  // namespace

}  // namespace peterson_r_graham_scan
