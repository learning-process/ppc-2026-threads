#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>

#include "peterson_r_graham_scan_seq/common/include/common.hpp"
#include "peterson_r_graham_scan_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace peterson_r_graham_scan_seq {

class PetersonRGrahamScanFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(PetersonRGrahamScanFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "circle_3"), std::make_tuple(5, "circle_5"),
                                            std::make_tuple(7, "circle_7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<PetersonRGrahamScanSeq, InType>(kTestParam, PPC_SETTINGS_peterson_r_graham_scan_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = PetersonRGrahamScanFuncTests::PrintFuncTestName<PetersonRGrahamScanFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, PetersonRGrahamScanFuncTests, kGtestValues, kPerfTestName);

void ExecutePipeline(const std::shared_ptr<PetersonRGrahamScanSeq> &task) {
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(PetersonRGrahamScanSeq, EmptyInput) {
  auto task = std::make_shared<PetersonRGrahamScanSeq>(0);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(PetersonRGrahamScanSeq, SinglePoint) {
  PointCloud pts = {Coordinate2D(5.0, 3.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonRGrahamScanSeq, TwoDistinctPoints) {
  PointCloud pts = {Coordinate2D(0.0, 0.0), Coordinate2D(3.0, 4.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(PetersonRGrahamScanSeq, CollinearPoints) {
  PointCloud pts = {Coordinate2D(0.0, 0.0), Coordinate2D(1.0, 0.0), Coordinate2D(2.0, 0.0), Coordinate2D(3.0, 0.0),
                    Coordinate2D(4.0, 0.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(PetersonRGrahamScanSeq, TrianglePoints) {
  PointCloud pts = {Coordinate2D(0.0, 0.0), Coordinate2D(4.0, 0.0), Coordinate2D(2.0, 3.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 3);
  EXPECT_EQ(static_cast<int>(task->GetHull().size()), 3);
}

TEST(PetersonRGrahamScanSeq, SquarePoints) {
  PointCloud pts = {Coordinate2D(0.0, 0.0), Coordinate2D(4.0, 0.0), Coordinate2D(4.0, 4.0), Coordinate2D(0.0, 4.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(PetersonRGrahamScanSeq, SquareWithInteriorPoint) {
  PointCloud pts = {Coordinate2D(0.0, 0.0), Coordinate2D(4.0, 0.0), Coordinate2D(4.0, 4.0), Coordinate2D(0.0, 4.0),
                    Coordinate2D(2.0, 2.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(PetersonRGrahamScanSeq, AllIdenticalPoints) {
  PointCloud pts = {Coordinate2D(3.0, 3.0), Coordinate2D(3.0, 3.0), Coordinate2D(3.0, 3.0), Coordinate2D(3.0, 3.0),
                    Coordinate2D(3.0, 3.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 1);
}

TEST(PetersonRGrahamScanSeq, NegativeInput) {
  auto task = std::make_shared<PetersonRGrahamScanSeq>(-1);
  EXPECT_FALSE(task->Validation());
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
}

TEST(PetersonRGrahamScanSeq, PointOnBoundary) {
  PointCloud pts = {Coordinate2D(0.0, 0.0), Coordinate2D(4.0, 0.0), Coordinate2D(2.0, 0.0), Coordinate2D(4.0, 4.0),
                    Coordinate2D(0.0, 4.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 4);
}

TEST(PetersonRGrahamScanSeq, VerticalCollinear) {
  PointCloud pts = {Coordinate2D(0.0, 0.0), Coordinate2D(0.0, 1.0), Coordinate2D(0.0, 2.0), Coordinate2D(0.0, 5.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 2);
}

TEST(PetersonRGrahamScanSeq, LargeCircle) {
  auto task = std::make_shared<PetersonRGrahamScanSeq>(1000);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 1000);
}

TEST(PetersonRGrahamScanSeq, HexagonWithCenter) {
  PointCloud pts = {Coordinate2D(2.0, 0.0),  Coordinate2D(1.0, 1.73),   Coordinate2D(-1.0, 1.73),
                    Coordinate2D(-2.0, 0.0), Coordinate2D(-1.0, -1.73), Coordinate2D(1.0, -1.73),
                    Coordinate2D(0.0, 0.0)};
  auto task = std::make_shared<PetersonRGrahamScanSeq>(static_cast<int>(pts.size()));
  task->SetPoints(pts);
  ExecutePipeline(task);
  EXPECT_EQ(task->GetOutput(), 6);
}

}  // namespace

}  // namespace peterson_r_graham_scan_seq
