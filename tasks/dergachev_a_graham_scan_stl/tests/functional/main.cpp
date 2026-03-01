#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "dergachev_a_graham_scan_stl/common/include/common.hpp"
#include "dergachev_a_graham_scan_stl/stl/include/ops_stl.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dergachev_a_graham_scan_stl {

class DergachevAGrahamScanFuncTestsStl : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (input_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(DergachevAGrahamScanFuncTestsStl, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "circle_3"), std::make_tuple(5, "circle_5"),
                                            std::make_tuple(7, "circle_7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<DergachevAGrahamScanSTL, InType>(kTestParam, PPC_SETTINGS_dergachev_a_graham_scan_stl));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = DergachevAGrahamScanFuncTestsStl::PrintFuncTestName<DergachevAGrahamScanFuncTestsStl>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, DergachevAGrahamScanFuncTestsStl, kGtestValues, kPerfTestName);

void RunstlPipeline(const std::shared_ptr<DergachevAGrahamScanSTL> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(DergachevAGrahamScanSTL, EmptyInput) {
  auto task = std::make_shared<DergachevAGrahamScanSTL>(0);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(DergachevAGrahamScanSTL, SinglePoint) {
  std::vector<Point> pts = {{.x = 5.0, .y = 3.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(DergachevAGrahamScanSTL, TwoDistinctPoints) {
  std::vector<Point> pts = {{.x = 0.0, .y = 0.0}, {.x = 3.0, .y = 4.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(DergachevAGrahamScanSTL, CollinearPoints) {
  std::vector<Point> pts = {
      {.x = 0.0, .y = 0.0}, {.x = 1.0, .y = 0.0}, {.x = 2.0, .y = 0.0}, {.x = 3.0, .y = 0.0}, {.x = 4.0, .y = 0.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(DergachevAGrahamScanSTL, TrianglePoints) {
  std::vector<Point> pts = {{.x = 0.0, .y = 0.0}, {.x = 4.0, .y = 0.0}, {.x = 2.0, .y = 3.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 3);
  EXPECT_EQ(static_cast<int>(task->GetHull().size()), 3);
}

TEST(DergachevAGrahamScanSTL, SquarePoints) {
  std::vector<Point> pts = {{.x = 0.0, .y = 0.0}, {.x = 4.0, .y = 0.0}, {.x = 4.0, .y = 4.0}, {.x = 0.0, .y = 4.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(DergachevAGrahamScanSTL, SquareWithInteriorPoint) {
  std::vector<Point> pts = {
      {.x = 0.0, .y = 0.0}, {.x = 4.0, .y = 0.0}, {.x = 4.0, .y = 4.0}, {.x = 0.0, .y = 4.0}, {.x = 2.0, .y = 2.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(DergachevAGrahamScanSTL, AllIdenticalPoints) {
  std::vector<Point> pts = {
      {.x = 3.0, .y = 3.0}, {.x = 3.0, .y = 3.0}, {.x = 3.0, .y = 3.0}, {.x = 3.0, .y = 3.0}, {.x = 3.0, .y = 3.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(DergachevAGrahamScanSTL, NegativeInput) {
  auto task = std::make_shared<DergachevAGrahamScanSTL>(-1);
  EXPECT_FALSE(task->Validation());
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(DergachevAGrahamScanSTL, PointOnBoundary) {
  std::vector<Point> pts = {
      {.x = 0.0, .y = 0.0}, {.x = 4.0, .y = 0.0}, {.x = 2.0, .y = 0.0}, {.x = 4.0, .y = 4.0}, {.x = 0.0, .y = 4.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(DergachevAGrahamScanSTL, VerticalCollinear) {
  std::vector<Point> pts = {{.x = 0.0, .y = 0.0}, {.x = 0.0, .y = 1.0}, {.x = 0.0, .y = 2.0}, {.x = 0.0, .y = 5.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(DergachevAGrahamScanSTL, LargeCircle) {
  auto task = std::make_shared<DergachevAGrahamScanSTL>(1000);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 1000);
}

TEST(DergachevAGrahamScanSTL, HexagonWithCenter) {
  std::vector<Point> pts = {{.x = 2.0, .y = 0.0},  {.x = 1.0, .y = 1.73},   {.x = -1.0, .y = 1.73},
                            {.x = -2.0, .y = 0.0}, {.x = -1.0, .y = -1.73}, {.x = 1.0, .y = -1.73},
                            {.x = 0.0, .y = 0.0}};
  auto task = std::make_shared<DergachevAGrahamScanSTL>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RunstlPipeline(task);
  EXPECT_EQ(task->GetOutput(), 6);
}

}  // namespace

}  // namespace dergachev_a_graham_scan_stl
