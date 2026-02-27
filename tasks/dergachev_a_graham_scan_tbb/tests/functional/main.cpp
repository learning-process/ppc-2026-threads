#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "dergachev_a_graham_scan_tbb/common/include/common.hpp"
#include "dergachev_a_graham_scan_tbb/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dergachev_a_graham_scan_tbb {

class DergachevAGrahamScanFuncTestsTbb : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(DergachevAGrahamScanFuncTestsTbb, MatmulFrtbbic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "circle_3"), std::make_tuple(5, "circle_5"),
                                            std::make_tuple(7, "circle_7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<DergachevAGrahamScanTBB, InType>(kTestParam, PPC_SETTINGS_dergachev_a_graham_scan_tbb));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = DergachevAGrahamScanFuncTestsTbb::PrintFuncTestName<DergachevAGrahamScanFuncTestsTbb>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, DergachevAGrahamScanFuncTestsTbb, kGtestValues, kPerfTestName);

void RuntbbPipeline(const std::shared_ptr<DergachevAGrahamScanTBB> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(DergachevAGrahamScanTBB, EmptyInput) {
  auto task = std::make_shared<DergachevAGrahamScanTBB>(0);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(DergachevAGrahamScanTBB, SinglePoint) {
  std::vector<Point> pts = {{.x = 5.0, .y = 3.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(DergachevAGrahamScanTBB, TwoDistinctPoints) {
  std::vector<Point> pts = {{.x = 0.0, .y = 0.0}, {.x = 3.0, .y = 4.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(DergachevAGrahamScanTBB, CollinearPoints) {
  std::vector<Point> pts = {
      {.x = 0.0, .y = 0.0}, {.x = 1.0, .y = 0.0}, {.x = 2.0, .y = 0.0}, {.x = 3.0, .y = 0.0}, {.x = 4.0, .y = 0.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(DergachevAGrahamScanTBB, TrianglePoints) {
  std::vector<Point> pts = {{.x = 0.0, .y = 0.0}, {.x = 4.0, .y = 0.0}, {.x = 2.0, .y = 3.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 3);
  EXPECT_EQ(static_cast<int>(task->GetHull().size()), 3);
}

TEST(DergachevAGrahamScanTBB, SquarePoints) {
  std::vector<Point> pts = {{.x = 0.0, .y = 0.0}, {.x = 4.0, .y = 0.0}, {.x = 4.0, .y = 4.0}, {.x = 0.0, .y = 4.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(DergachevAGrahamScanTBB, SquareWithInteriorPoint) {
  std::vector<Point> pts = {
      {.x = 0.0, .y = 0.0}, {.x = 4.0, .y = 0.0}, {.x = 4.0, .y = 4.0}, {.x = 0.0, .y = 4.0}, {.x = 2.0, .y = 2.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(DergachevAGrahamScanTBB, AllIdenticalPoints) {
  std::vector<Point> pts = {
      {.x = 3.0, .y = 3.0}, {.x = 3.0, .y = 3.0}, {.x = 3.0, .y = 3.0}, {.x = 3.0, .y = 3.0}, {.x = 3.0, .y = 3.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(DergachevAGrahamScanTBB, NegativeInput) {
  auto task = std::make_shared<DergachevAGrahamScanTBB>(-1);
  EXPECT_FALSE(task->Validation());
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(DergachevAGrahamScanTBB, PointOnBoundary) {
  std::vector<Point> pts = {
      {.x = 0.0, .y = 0.0}, {.x = 4.0, .y = 0.0}, {.x = 2.0, .y = 0.0}, {.x = 4.0, .y = 4.0}, {.x = 0.0, .y = 4.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(DergachevAGrahamScanTBB, VerticalCollinear) {
  std::vector<Point> pts = {{.x = 0.0, .y = 0.0}, {.x = 0.0, .y = 1.0}, {.x = 0.0, .y = 2.0}, {.x = 0.0, .y = 5.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(DergachevAGrahamScanTBB, LargeCircle) {
  auto task = std::make_shared<DergachevAGrahamScanTBB>(1000);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 1000);
}

TEST(DergachevAGrahamScanTBB, HexagonWithCenter) {
  std::vector<Point> pts = {{.x = 2.0, .y = 0.0},  {.x = 1.0, .y = 1.73},   {.x = -1.0, .y = 1.73},
                            {.x = -2.0, .y = 0.0}, {.x = -1.0, .y = -1.73}, {.x = 1.0, .y = -1.73},
                            {.x = 0.0, .y = 0.0}};
  auto task = std::make_shared<DergachevAGrahamScanTBB>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  RuntbbPipeline(task);
  EXPECT_EQ(task->GetOutput(), 6);
}

}  // namespace

}  // namespace dergachev_a_graham_scan_tbb
