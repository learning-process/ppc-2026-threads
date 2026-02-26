#include <gtest/gtest.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <vector>

#include "kondrashova_v_marking_components_seq/common/include/common.hpp"
#include "kondrashova_v_marking_components_seq/seq/include/ops_seq.hpp"
#include "performance/include/performance.hpp"
#include "util/include/perf_test_util.hpp"

namespace kondrashova_v_marking_components_seq {

namespace {
void SetTimer(ppc::performance::PerfAttr &perf_attrs) {
  const auto t0 = std::chrono::high_resolution_clock::now();
  perf_attrs.current_timer = [t0] {
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = now - t0;
    return diff.count();
  };
}

bool CheckLabels(const OutType &output_data, int expected_count, int size) {
  if (output_data.count != expected_count) return false;
  if (output_data.labels.size() != size) return false;
  if (!output_data.labels.empty() && (output_data.labels[0].size() != size))
    return false;
  return true;
}
}  // namespace

class AllOnesPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    SetTimer(perf_attrs);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    return CheckLabels(output_data, 0, kSize);
  }
  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(static_cast<size_t>(kSize * kSize), 1);
    return input;
  }
};

class AllZerosPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    SetTimer(perf_attrs);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    return CheckLabels(output_data, 1, kSize);
  }
  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(static_cast<size_t>(kSize * kSize), 0);
    return input;
  }
};

class ChessboardPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    SetTimer(perf_attrs);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count <= 0) return false;
    if (output_data.labels.size() != kSize) return false;
    if (!output_data.labels.empty() && (output_data.labels[0].size() != kSize))
      return false;
    return true;
  }
  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.resize(static_cast<size_t>(kSize * kSize));
    for (int row = 0; row < kSize; ++row) {
      for (int col = 0; col < kSize; ++col) {
        input.data[static_cast<size_t>((row * kSize) + col)] =
            static_cast<uint8_t>((row + col) % 2);
      }
    }
    return input;
  }
};

class SparseDotsPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    SetTimer(perf_attrs);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count <= 0) return false;
    if (output_data.labels.size() != kSize) return false;
    if (!output_data.labels.empty() && (output_data.labels[0].size() != kSize))
      return false;
    return true;
  }
  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(static_cast<size_t>(kSize * kSize), 1);
    for (int row = 0; row < kSize; row += 8) {
      for (int col = 0; col < kSize; col += 8) {
        input.data[static_cast<size_t>((row * kSize) + col)] = 0;
      }
    }
    return input;
  }
};

class StripesPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    SetTimer(perf_attrs);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count <= 0) return false;
    if (output_data.labels.size() != kSize) return false;
    if (!output_data.labels.empty() && (output_data.labels[0].size() != kSize))
      return false;
    return true;
  }
  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(static_cast<size_t>(kSize * kSize), 1);
    for (int row = 0; row < kSize; row += 4) {
      for (int col = 0; col < kSize; ++col) {
        input.data[static_cast<size_t>((row * kSize) + col)] = 0;
      }
    }
    return input;
  }
};

class BlocksPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;
  static constexpr int kBlock = 32;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    SetTimer(perf_attrs);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count <= 0) return false;
    if (output_data.labels.size() != kSize) return false;
    if (!output_data.labels.empty() && (output_data.labels[0].size() != kSize))
      return false;
    return true;
  }
  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(static_cast<size_t>(kSize * kSize), 1);
    for (int by = 0; by < kSize; by += kBlock * 2) {
      for (int bx = 0; bx < kSize; bx += kBlock * 2) {
        for (int row = 0; row < kBlock && by + row < kSize; ++row) {
          for (int col = 0; col < kBlock && bx + col < kSize; ++col) {
            input.data[static_cast<size_t>(((by + row) * kSize) + (bx + col))] = 0;
          }
        }
      }
    }
    return input;
  }
};

TEST_P(AllOnesPerfTest, RunPerfModes) { ExecuteTest(GetParam()); }
TEST_P(AllZerosPerfTest, RunPerfModes) { ExecuteTest(GetParam()); }
TEST_P(ChessboardPerfTest, RunPerfModes) { ExecuteTest(GetParam()); }
TEST_P(SparseDotsPerfTest, RunPerfModes) { ExecuteTest(GetParam()); }
TEST_P(StripesPerfTest, RunPerfModes) { ExecuteTest(GetParam()); }
TEST_P(BlocksPerfTest, RunPerfModes) { ExecuteTest(GetParam()); }

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KondrashovaVTaskSEQ>(
        PPC_SETTINGS_kondrashova_v_marking_components_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

INSTANTIATE_TEST_SUITE_P(AllOnes_RunModeTests, AllOnesPerfTest, kGtestValues,
                         AllOnesPerfTest::CustomPerfTestName);
INSTANTIATE_TEST_SUITE_P(AllZeros_RunModeTests, AllZerosPerfTest, kGtestValues,
                         AllZerosPerfTest::CustomPerfTestName);
INSTANTIATE_TEST_SUITE_P(Chessboard_RunModeTests, ChessboardPerfTest,
                         kGtestValues, ChessboardPerfTest::CustomPerfTestName);
INSTANTIATE_TEST_SUITE_P(SparseDots_RunModeTests, SparseDotsPerfTest,
                         kGtestValues, SparseDotsPerfTest::CustomPerfTestName);
INSTANTIATE_TEST_SUITE_P(Stripes_RunModeTests, StripesPerfTest, kGtestValues,
                         StripesPerfTest::CustomPerfTestName);
INSTANTIATE_TEST_SUITE_P(Blocks_RunModeTests, BlocksPerfTest, kGtestValues,
                         BlocksPerfTest::CustomPerfTestName);

}  // namespace

}  // namespace kondrashova_v_marking_components_seq