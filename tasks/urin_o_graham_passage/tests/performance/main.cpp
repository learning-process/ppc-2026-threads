#include <gtest/gtest.h>

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

namespace urin_o_graham_passage {

static bool IsConvexHull(const std::vector<Point> &hull) {
  if (hull.size() < 3) {
    return true;
  }

  for (size_t i = 0; i < hull.size(); i++) {
    size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
    size_t next = (i + 1) % hull.size();

    if (UrinOGrahamPassageSEQ::Orientation(hull[prev], hull[i], hull[next]) < 0) {
      return false;
    }
  }
  return true;
}

class UrinOGrahamPassagePerfTest : public ::testing::Test {
 protected:
  InType GenerateRandomPoints(size_t num_points) {
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

TEST_F(UrinOGrahamPassagePerfTest, SeqPerformance) {
  const size_t num_points = 10000;
  InType input_points = GenerateRandomPoints(num_points);

  UrinOGrahamPassageSEQ task(input_points);

  auto start = std::chrono::high_resolution_clock::now();

  EXPECT_TRUE(task.Validation());
  EXPECT_TRUE(task.PreProcessing());
  EXPECT_TRUE(task.Run());
  EXPECT_TRUE(task.PostProcessing());

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const auto &hull = task.GetOutput();
  // ИСПРАВЛЕНО: сравниваем size_t с size_t
  EXPECT_GE(hull.size(), size_t(3));
  EXPECT_TRUE(IsConvexHull(hull));

  std::cout << "SEQ version with " << num_points << " points took " << duration.count() << " ms" << std::endl;
  std::cout << "Convex hull size: " << hull.size() << std::endl;
}

TEST_F(UrinOGrahamPassagePerfTest, DifferentSizes) {
  std::vector<size_t> sizes = {100, 500, 1000, 5000, 10000};

  std::cout << "\nPerformance test with different sizes:" << std::endl;

  for (size_t size : sizes) {
    InType test_points = GenerateRandomPoints(size);

    UrinOGrahamPassageSEQ task(test_points);

    auto start = std::chrono::high_resolution_clock::now();

    task.Validation();
    task.PreProcessing();
    task.Run();
    task.PostProcessing();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    const auto &hull = task.GetOutput();
    // ИСПРАВЛЕНО: все сравнения с size_t
    if (hull.size() >= size_t(3)) {
      std::cout << "Size " << size << ": " << duration.count() << " ms, "
                << "hull size: " << hull.size() << std::endl;
    } else {
      std::cout << "Size " << size << ": " << duration.count() << " ms, "
                << "hull size: " << hull.size() << " (invalid)" << std::endl;
    }
  }
}

}  // namespace urin_o_graham_passage
