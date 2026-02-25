#include <gtest/gtest.h>

#include <chrono>
#include <vector>

#include "kondrashova_v_marking_components_seq/common/include/common.hpp"
#include "kondrashova_v_marking_components_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kondrashova_v_marking_components_seq {

// 1. Лучший случай — всё единицы (фон, 0 компонент)
class AllOnesPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto t0 = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [t0] {  // захват по ЗНАЧЕНИЮ!
      auto now = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = now - t0;
      return diff.count();
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count != 0) {
      return false;
    }
    if (static_cast<int>(output_data.labels.size()) != kSize) {
      return false;
    }
    if (!output_data.labels.empty() && static_cast<int>(output_data.labels[0].size()) != kSize) {
      return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(kSize * kSize, 1);
    return input;
  }
};

// 2. Всё нули — одна большая компонента
class AllZerosPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto t0 = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [t0] {
      auto now = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = now - t0;
      return diff.count();
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count != 1) {
      return false;
    }
    if (static_cast<int>(output_data.labels.size()) != kSize) {
      return false;
    }
    if (!output_data.labels.empty() && static_cast<int>(output_data.labels[0].size()) != kSize) {
      return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(kSize * kSize, 0);
    return input;
  }
};

// 3. Шахматная доска (максимум компонент)
class ChessboardPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto t0 = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [t0] {
      auto now = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = now - t0;
      return diff.count();
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count <= 0) {
      return false;
    }
    if (static_cast<int>(output_data.labels.size()) != kSize) {
      return false;
    }
    if (!output_data.labels.empty() && static_cast<int>(output_data.labels[0].size()) != kSize) {
      return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.resize(kSize * kSize);

    for (int y = 0; y < kSize; ++y) {
      for (int x = 0; x < kSize; ++x) {
        input.data[y * kSize + x] = static_cast<uint8_t>((x + y) % 2);
      }
    }
    return input;
  }
};

// 4. Разреженные точки
class SparseDotsPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto t0 = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [t0] {
      auto now = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = now - t0;
      return diff.count();
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count <= 0) {
      return false;
    }
    if (static_cast<int>(output_data.labels.size()) != kSize) {
      return false;
    }
    if (!output_data.labels.empty() && static_cast<int>(output_data.labels[0].size()) != kSize) {
      return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(kSize * kSize, 1);

    for (int y = 0; y < kSize; y += 8) {
      for (int x = 0; x < kSize; x += 8) {
        input.data[y * kSize + x] = 0;
      }
    }
    return input;
  }
};

// 5. Полосы
class StripesPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto t0 = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [t0] {
      auto now = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = now - t0;
      return diff.count();
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count <= 0) {
      return false;
    }
    if (static_cast<int>(output_data.labels.size()) != kSize) {
      return false;
    }
    if (!output_data.labels.empty() && static_cast<int>(output_data.labels[0].size()) != kSize) {
      return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(kSize * kSize, 1);

    for (int y = 0; y < kSize; y += 4) {
      for (int x = 0; x < kSize; ++x) {
        input.data[y * kSize + x] = 0;
      }
    }
    return input;
  }
};

// 6. Крупные блоки
class BlocksPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 512;
  static constexpr int kBlock = 32;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto t0 = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [t0] {
      auto now = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = now - t0;
      return diff.count();
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count <= 0) {
      return false;
    }
    if (static_cast<int>(output_data.labels.size()) != kSize) {
      return false;
    }
    if (!output_data.labels.empty() && static_cast<int>(output_data.labels[0].size()) != kSize) {
      return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    InType input;
    input.width = kSize;
    input.height = kSize;
    input.data.assign(kSize * kSize, 1);

    for (int by = 0; by < kSize; by += kBlock * 2) {
      for (int bx = 0; bx < kSize; bx += kBlock * 2) {
        for (int y = 0; y < kBlock && by + y < kSize; ++y) {
          for (int x = 0; x < kBlock && bx + x < kSize; ++x) {
            input.data[(by + y) * kSize + (bx + x)] = 0;
          }
        }
      }
    }
    return input;
  }
};

// ---- Запуск тестов ----

TEST_P(AllOnesPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}
TEST_P(AllZerosPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}
TEST_P(ChessboardPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}
TEST_P(SparseDotsPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}
TEST_P(StripesPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}
TEST_P(BlocksPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KondrashovaVTaskSEQ>(PPC_SETTINGS_kondrashova_v_marking_components_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

INSTANTIATE_TEST_SUITE_P(AllOnes_RunModeTests, AllOnesPerfTest, kGtestValues, AllOnesPerfTest::CustomPerfTestName);

INSTANTIATE_TEST_SUITE_P(AllZeros_RunModeTests, AllZerosPerfTest, kGtestValues, AllZerosPerfTest::CustomPerfTestName);

INSTANTIATE_TEST_SUITE_P(Chessboard_RunModeTests, ChessboardPerfTest, kGtestValues,
                         ChessboardPerfTest::CustomPerfTestName);

INSTANTIATE_TEST_SUITE_P(SparseDots_RunModeTests, SparseDotsPerfTest, kGtestValues,
                         SparseDotsPerfTest::CustomPerfTestName);

INSTANTIATE_TEST_SUITE_P(Stripes_RunModeTests, StripesPerfTest, kGtestValues, StripesPerfTest::CustomPerfTestName);

INSTANTIATE_TEST_SUITE_P(Blocks_RunModeTests, BlocksPerfTest, kGtestValues, BlocksPerfTest::CustomPerfTestName);

}  // namespace

}  // namespace kondrashova_v_marking_components_seq
