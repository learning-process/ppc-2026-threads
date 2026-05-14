#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <tuple>
#include <utility>

#include "klimovich_v_crs_complex_mat_mul/common/include/common.hpp"
#include "klimovich_v_crs_complex_mat_mul/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace klimovich_v_crs_complex_mat_mul {
namespace {

CrsMatrix BuildScatteredMatrix(int n, std::uint32_t step, std::uint32_t shift) {
  CrsMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    const auto u = static_cast<std::uint32_t>(i);
    const int c1 = static_cast<int>((u * step + shift) % static_cast<std::uint32_t>(n));
    int c2 = static_cast<int>((u * (step + 13U) + shift + 7U) % static_cast<std::uint32_t>(n));
    if (c2 == c1) {
      c2 = (c2 + 1) % n;
    }
    int low = c1;
    int high = c2;
    if (low > high) {
      std::swap(low, high);
    }
    m.col_indices.push_back(low);
    m.data.emplace_back(static_cast<double>((i % 7) + 1), static_cast<double>((i % 5) - 2));
    m.col_indices.push_back(high);
    m.data.emplace_back(static_cast<double>((i % 11) - 5), static_cast<double>((i % 9) + 1));
    m.row_offsets[i + 1] = static_cast<int>(m.col_indices.size());
  }
  return m;
}

}  // namespace

class KlimovichVCrsComplexPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSize = 1500;
  InType input_data_;

 public:
  void SetUp() override {
    auto lhs = BuildScatteredMatrix(kSize, 17U, 3U);
    auto rhs = BuildScatteredMatrix(kSize, 23U, 11U);
    input_data_ = std::make_tuple(lhs, rhs);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    if (ppc::util::IsUnderMpirun()) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    if (rank != 0) {
      return true;
    }
    return output_data.n_rows == kSize && output_data.n_cols == kSize;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KlimovichVCrsComplexPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KlimovichVCrsComplexMatMulSeq>(PPC_SETTINGS_klimovich_v_crs_complex_mat_mul);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KlimovichVCrsComplexPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KlimovichVCrsComplexPerfTest, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace klimovich_v_crs_complex_mat_mul
