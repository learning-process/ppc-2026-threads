#include <gtest/gtest.h>

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

namespace urin_o_graham_passage {

class UrinOGrahamPassagePerfTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Генерируем точки в SetUp, чтобы можно было использовать в разных тестах
  }

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

  bool IsConvexHull(const OutType &hull) {
    if (hull.size() < 3) {
      return false;
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
};

TEST_F(UrinOGrahamPassagePerfTest, SeqPerformance) {
  const size_t num_points = 10000;
  InType input_points = GenerateRandomPoints(num_points);

  UrinOGrahamPassageSEQ task(input_points);
  ppc::task::Task<InType, OutType> *base_task = &task;

  auto start = std::chrono::high_resolution_clock::now();

  EXPECT_TRUE(base_task->Validation());
  EXPECT_TRUE(base_task->PreProcessing());
  EXPECT_TRUE(base_task->Run());
  EXPECT_TRUE(base_task->PostProcessing());

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const auto &hull = task.GetOutput();
  EXPECT_TRUE(IsConvexHull(hull));
  EXPECT_GE(hull.size(), 3);

  std::cout << "SEQ version with " << num_points << " points took " << duration.count() << " ms" << std::endl;
  std::cout << "Convex hull size: " << hull.size() << std::endl;
}

TEST_F(UrinOGrahamPassagePerfTest, DifferentSizes) {
  std::vector<size_t> sizes = {100, 500, 1000, 5000, 10000};

  std::cout << "\nPerformance test with different sizes:" << std::endl;

  for (size_t size : sizes) {
    InType test_points = GenerateRandomPoints(size);

    UrinOGrahamPassageSEQ task(test_points);
    ppc::task::Task<InType, OutType> *base_task = &task;

    auto start = std::chrono::high_resolution_clock::now();

    base_task->Validation();
    base_task->PreProcessing();
    base_task->Run();
    base_task->PostProcessing();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Size " << size << ": " << duration.count() << " ms, "
              << "hull size: " << task.GetOutput().size() << std::endl;
  }
}

}  // namespace urin_o_graham_passage
