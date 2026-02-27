#include <gtest/gtest.h>
#include <random>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"
#include "kutuzov_i_convex_hull_jarvis/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kutuzov_i_convex_hull_jarvis {

class KutuzovIRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 5000;
  InType input_data_ = {};

  void SetUp() override {
    input_data_ = {};

    std::mt19937 rng(1);
    std::uniform_real_distribution<double> dist(-1000, 1000);

    for (int i = 0; i < kCount_; i++) {
      double random_x = dist(rng);
      double random_y = dist(rng);

      input_data_.push_back({random_x, random_y});
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    OutType expected_output = {};

    if (input_data_.size() < 3) {
      expected_output = input_data_;
      return true;
    }
    
    // Finding left-most point
    int leftmost = 0;
    double leftmost_x = std::get<0>(input_data_[leftmost]);
    double leftmost_y = std::get<1>(input_data_[leftmost]);

    for (int i = 0; i <  static_cast<int>(input_data_.size()); i++) {
      double x = std::get<0>(input_data_[i]);
      double y = std::get<1>(input_data_[i]);

      if ((x < leftmost_x) || ((x == leftmost_x) && (y < leftmost_y))) {
        leftmost = i;
        leftmost_x = std::get<0>(input_data_[leftmost]);
        leftmost_y = std::get<1>(input_data_[leftmost]);
      }
    }

    // Main loop
    int current = leftmost;
    double current_x = std::get<0>(input_data_[current]);
    double current_y = std::get<1>(input_data_[current]);
    
    const double epsilon = 1e-9;

    do {
      // Adding current point to the hull 
      expected_output.push_back(input_data_[current]);

      // Finding the next point of the hull
      int next = (current + 1) %  static_cast<int>(input_data_.size());
      double next_x = std::get<0>(input_data_[next]);
      double next_y = std::get<1>(input_data_[next]);
      
      for (int i = 0; i < static_cast<int>(input_data_.size()); i++) {
        if (i == current) continue;
        
        double i_x = std::get<0>(input_data_[i]);
        double i_y = std::get<1>(input_data_[i]);

        double cross = CrossProduct(current_x, current_y, next_x, next_y, i_x, i_y);

        if (cross < -epsilon || ((std::abs(cross) < epsilon) && (DistanceSquared(current_x, current_y, i_x, i_y) > DistanceSquared(current_x, current_y, next_x, next_y)))) {
          next = i;
          next_x = std::get<0>(input_data_[next]);
          next_y = std::get<1>(input_data_[next]);
        }
      }

      current = next;
      current_x = next_x;
      current_y = next_y;

    } while (current != leftmost); // Loop until we wrap around to the first point

    return expected_output == output_data;
  }

  double DistanceSquared(double a_x, double a_y, double b_x, double b_y) {
    return (a_x - b_x) * (a_x - b_x) + (a_y - b_y) * (a_y - b_y);
  }

  double CrossProduct(double o_x, double o_y, double a_x, double a_y, double b_x, double b_y) {
    return (a_x - o_x) * (b_y - o_y) - (a_y - o_y) * (b_x - o_x);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KutuzovIRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KutuzovITestConvexHullSEQ>(PPC_SETTINGS_kutuzov_i_convex_hull_jarvis);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KutuzovIRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KutuzovIRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kutuzov_i_convex_hull_jarvis
