#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/omp/include/ops_omp.hpp"  // ИЗМЕНЕНО

namespace urin_o_graham_passage {
namespace {

bool IsConvexHull(const std::vector<Point> &hull) {
  if (hull.size() < 3) {
    return true;
  }

  for (size_t i = 0; i < hull.size(); ++i) {
    size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
    size_t next = (i + 1) % hull.size();

    if (UrinOGrahamPassageOMP::Orientation(hull[prev], hull[i], hull[next]) < 0) {  // ИЗМЕНЕНО
      return false;
    }
  }
  return true;
}

class UrinOGrahamPassagePerfTest : public ::testing::Test {
 protected:
  static InType GenerateRandomPoints(size_t num_points) {
    InType points;
    points.reserve(num_points);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);

    for (size_t i = 0; i < num_points; ++i) {
      points.emplace_back(dist(gen), dist(gen));
    }

    return points;
  }
};

bool ValidateTask(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  return task.Validation();
}

bool PreProcessTask(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  return task.PreProcessing();
}

bool RunTask(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  return task.Run();
}

bool PostProcessTask(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  return task.PostProcessing();
}

void ExpectValidation(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  EXPECT_TRUE(ValidateTask(task));
}

void ExpectPreProcessing(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  EXPECT_TRUE(PreProcessTask(task));
}

void ExpectRun(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  EXPECT_TRUE(RunTask(task));
}

void ExpectPostProcessing(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  EXPECT_TRUE(PostProcessTask(task));
}

void RunTaskPipeline(UrinOGrahamPassageOMP &task) {  // ИЗМЕНЕНО
  ExpectValidation(task);
  ExpectPreProcessing(task);
  ExpectRun(task);
  ExpectPostProcessing(task);
}

void CheckHullValidity(const std::vector<Point> &hull) {
  EXPECT_GE(hull.size(), static_cast<size_t>(3));
  EXPECT_TRUE(IsConvexHull(hull));
}

void PrintPerformanceResult(size_t num_points, int64_t ms, size_t hull_size) {
  std::cout << "OMP version with " << num_points << " points took " << ms << " ms\n";  // ИЗМЕНЕНО
  std::cout << "Convex hull size: " << hull_size << "\n";
}

TEST_F(UrinOGrahamPassagePerfTest, OmpPerformance) {  // ИЗМЕНЕНО
  const size_t num_points = 10000;
  InType input_points = GenerateRandomPoints(num_points);

  UrinOGrahamPassageOMP task(input_points);  // ИЗМЕНЕНО

  auto start = std::chrono::high_resolution_clock::now();
  RunTaskPipeline(task);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const auto &hull = task.GetOutput();
  CheckHullValidity(hull);
  PrintPerformanceResult(num_points, static_cast<int64_t>(duration.count()), hull.size());
}

TEST_F(UrinOGrahamPassagePerfTest, DifferentSizes) {
  std::vector<size_t> sizes = {100, 500, 1000, 5000, 10000};

  std::cout << "\nPerformance test with different sizes:\n";

  for (size_t size : sizes) {
    InType test_points = GenerateRandomPoints(size);
    UrinOGrahamPassageOMP task(test_points);  // ИЗМЕНЕНО

    auto start = std::chrono::high_resolution_clock::now();
    RunTaskPipeline(task);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    const auto &hull = task.GetOutput();

    if (hull.size() >= static_cast<size_t>(3)) {
      std::cout << "Size " << size << ": " << duration.count() << " ms, "
                << "hull size: " << hull.size() << "\n";
    } else {
      std::cout << "Size " << size << ": " << duration.count() << " ms, "
                << "hull size: " << hull.size() << " (invalid)\n";
    }
  }
}

}  // namespace
}  // namespace urin_o_graham_passage
