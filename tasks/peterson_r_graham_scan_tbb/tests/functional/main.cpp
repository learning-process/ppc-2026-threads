#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "peterson_r_graham_scan_tbb/common/include/common.hpp"
#include "peterson_r_graham_scan_tbb/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace peterson_r_graham_scan_tbb {

class PetersonGrahamScannerTBBFuncTests : public ppc::util::BaseRunFuncTests<InputValue, OutputValue, TestParameters> {
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

TEST_P(PetersonGrahamScannerTBBFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestParameters, 3> kTestCases = {std::make_tuple(3, "circle_3"), std::make_tuple(5, "circle_5"),
                                                  std::make_tuple(7, "circle_7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<PetersonGrahamScannerTBB, InputValue>(kTestCases, PPC_SETTINGS_peterson_r_graham_scan_tbb));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestNameGenerator = PetersonGrahamScannerTBBFuncTests::PrintFuncTestName<PetersonGrahamScannerTBBFuncTests>;

INSTANTIATE_TEST_SUITE_P(DefaultTests, PetersonGrahamScannerTBBFuncTests, kGtestValues, kTestNameGenerator);

void ExecutePipeline(const std::shared_ptr<PetersonGrahamScannerTBB> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(PetersonGrahamScannerTBB, EmptyInput) {
  auto task = std::make_shared<PetersonGrahamScannerTBB>(0);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(PetersonGrahamScannerTBB, SinglePoint) {
  std::vector<Point> pts = {Point(5.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerTBB, TwoDistinctPoints) {
  std::vector<Point> pts = {Point(0.0, 0.0), Point(3.0, 4.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(PetersonGrahamScannerTBB, CollinearPoints) {
  std::vector<Point> pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0), Point(4.0, 0.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(PetersonGrahamScannerTBB, TrianglePoints) {
  std::vector<Point> pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 3);
  EXPECT_EQ(static_cast<int>(task->GetConvexHull().size()), 3);
}

TEST(PetersonGrahamScannerTBB, SquarePoints) {
  std::vector<Point> pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(PetersonGrahamScannerTBB, SquareWithInteriorPoint) {
  std::vector<Point> pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0), Point(2.0, 2.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(PetersonGrahamScannerTBB, AllIdenticalPoints) {
  std::vector<Point> pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonGrahamScannerTBB, NegativeInput) {
  auto task = std::make_shared<PetersonGrahamScannerTBB>(-1);
  EXPECT_FALSE(task->Validation());
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(PetersonGrahamScannerTBB, PointOnBoundary) {
  std::vector<Point> pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(PetersonGrahamScannerTBB, VerticalCollinear) {
  std::vector<Point> pts = {Point(0.0, 0.0), Point(0.0, 1.0), Point(0.0, 2.0), Point(0.0, 5.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(PetersonGrahamScannerTBB, LargeCircle) {
  auto task = std::make_shared<PetersonGrahamScannerTBB>(1000);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 1000);
}

TEST(PetersonGrahamScannerTBB, HexagonWithCenter) {
  std::vector<Point> pts = {Point(2.0, 0.0),    Point(1.0, 1.73),  Point(-1.0, 1.73), Point(-2.0, 0.0),
                            Point(-1.0, -1.73), Point(1.0, -1.73), Point(0.0, 0.0)};
  auto task = std::make_shared<PetersonGrahamScannerTBB>(static_cast<int>(pts.size()));
  task->LoadPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 6);
}

}  // namespace

}  // namespace peterson_r_graham_scan_tbb
