#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string_view>
#include <vector>

#include "urin_o_graham_passage/all/include/ops_all.hpp"
#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/omp/include/ops_omp.hpp"
#include "urin_o_graham_passage/seq/include/ops_seq.hpp"
#include "urin_o_graham_passage/stl/include/ops_stl.hpp"
#include "urin_o_graham_passage/tbb/include/ops_tbb.hpp"

namespace urin_o_graham_passage {
namespace {

double Orientation(const Point &p, const Point &q, const Point &r) {
  return ((q.x - p.x) * (r.y - p.y)) - ((q.y - p.y) * (r.x - p.x));
}

bool IsConvexHull(const std::vector<Point> &hull) {
  if (hull.size() < 3) {
    return true;
  }

  for (size_t i = 0; i < hull.size(); ++i) {
    size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
    size_t next = (i + 1) % hull.size();

    if (Orientation(hull[prev], hull[i], hull[next]) < -1e-10) {
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

template <class TaskType>
bool ValidateTask(TaskType &task) {
  return task.Validation();
}

template <class TaskType>
bool PreProcessTask(TaskType &task) {
  return task.PreProcessing();
}

template <class TaskType>
bool RunTask(TaskType &task) {
  return task.Run();
}

template <class TaskType>
bool PostProcessTask(TaskType &task) {
  return task.PostProcessing();
}

template <class TaskType>
void RunTaskPipeline(TaskType &task) {
  EXPECT_TRUE(ValidateTask(task));
  EXPECT_TRUE(PreProcessTask(task));
  EXPECT_TRUE(RunTask(task));
  EXPECT_TRUE(PostProcessTask(task));
}

void CheckHullValidity(const std::vector<Point> &hull) {
  EXPECT_GE(hull.size(), static_cast<size_t>(3));
  EXPECT_TRUE(IsConvexHull(hull));
}

void PrintPerformanceResult(std::string_view version, size_t num_points, int64_t ms, size_t hull_size) {
  std::cout << version << " version with " << num_points << " points took " << ms << " ms\n";
  std::cout << "Convex hull size: " << hull_size << "\n";
}

template <class TaskType>
void RunPerformanceCase(std::string_view version, const InType &input_points) {
  TaskType task(input_points);

  auto start = std::chrono::high_resolution_clock::now();
  RunTaskPipeline(task);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const auto &hull = task.GetOutput();
  CheckHullValidity(hull);
  PrintPerformanceResult(version, input_points.size(), static_cast<int64_t>(duration.count()), hull.size());
}

template <class TaskType>
void RunDifferentSizeCase(std::string_view version, const InType &test_points) {
  TaskType task(test_points);

  auto start = std::chrono::high_resolution_clock::now();
  RunTaskPipeline(task);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const auto &hull = task.GetOutput();

  if (hull.size() >= static_cast<size_t>(3)) {
    std::cout << version << " size " << test_points.size() << ": " << duration.count() << " ms, "
              << "hull size: " << hull.size() << "\n";
  } else {
    std::cout << version << " size " << test_points.size() << ": " << duration.count() << " ms, "
              << "hull size: " << hull.size() << " (invalid)\n";
  }
}

TEST_F(UrinOGrahamPassagePerfTest, Performance) {
  const size_t num_points = 10000;
  InType input_points = GenerateRandomPoints(num_points);

  RunPerformanceCase<UrinOGrahamPassageSEQ>("SEQ", input_points);
  RunPerformanceCase<UrinOGrahamPassageSTL>("STL", input_points);
  RunPerformanceCase<UrinOGrahamPassageOMP>("OMP", input_points);
  RunPerformanceCase<UrinOGrahamPassageTBB>("TBB", input_points);
  RunPerformanceCase<UrinOGrahamPassageALL>("ALL", input_points);
}

TEST_F(UrinOGrahamPassagePerfTest, DifferentSizes) {
  std::vector<size_t> sizes = {100, 500, 1000, 5000, 10000};

  std::cout << "\nPerformance test with different sizes:\n";

  for (size_t size : sizes) {
    InType test_points = GenerateRandomPoints(size);
    RunDifferentSizeCase<UrinOGrahamPassageSEQ>("SEQ", test_points);
    RunDifferentSizeCase<UrinOGrahamPassageSTL>("STL", test_points);
    RunDifferentSizeCase<UrinOGrahamPassageOMP>("OMP", test_points);
    RunDifferentSizeCase<UrinOGrahamPassageTBB>("TBB", test_points);
    RunDifferentSizeCase<UrinOGrahamPassageALL>("ALL", test_points);
  }
}

}  // namespace
}  // namespace urin_o_graham_passage
