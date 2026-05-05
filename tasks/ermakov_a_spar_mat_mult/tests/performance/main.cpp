#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "ermakov_a_spar_mat_mult/all/include/ops_all.hpp"
#include "ermakov_a_spar_mat_mult/common/include/common.hpp"
#include "ermakov_a_spar_mat_mult/omp/include/ops_omp.hpp"
#include "ermakov_a_spar_mat_mult/seq/include/ops_seq.hpp"
#include "ermakov_a_spar_mat_mult/stl/include/ops_stl.hpp"
#include "ermakov_a_spar_mat_mult/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace ermakov_a_spar_mat_mult {

namespace {

constexpr std::uint32_t kPerfSeedA = 0x13579BDFU;
constexpr std::uint32_t kPerfSeedB = 0x2468ACE0U;

MatrixCRS MakeRandomCRS(int n, double density, std::uint32_t seed) {
  MatrixCRS matrix;
  matrix.rows = n;
  matrix.cols = n;
  matrix.row_ptr.resize(static_cast<std::size_t>(n) + 1ULL, 0);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis_val(-5.0, 5.0);
  std::uniform_real_distribution<double> dis_prob(0.0, 1.0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (dis_prob(gen) < density) {
        matrix.values.emplace_back(dis_val(gen), dis_val(gen));
        matrix.col_index.push_back(j);
      }
    }
    matrix.row_ptr[static_cast<std::size_t>(i) + 1ULL] = static_cast<int>(matrix.values.size());
  }
  return matrix;
}

}  // namespace

class ErmakovARunPerfTestSparMatMult : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int k_count = 15000;
  InType input_data{};

  void SetUp() override {
    const double density = 0.001;

    input_data.A = MakeRandomCRS(k_count, density, kPerfSeedA);
    input_data.B = MakeRandomCRS(k_count, density, kPerfSeedB);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != k_count || output_data.cols != k_count) {
      return false;
    }

    if (output_data.row_ptr.size() != static_cast<size_t>(output_data.rows) + 1) {
      return false;
    }

    if (output_data.values.size() != output_data.col_index.size()) {
      return false;
    }

    if (output_data.row_ptr.front() != 0) {
      return false;
    }

    if (static_cast<size_t>(output_data.row_ptr.back()) != output_data.values.size()) {
      return false;
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(ErmakovARunPerfTestSparMatMult, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ErmakovASparMatMultSEQ, ErmakovASparMatMultOMP, ErmakovASparMatMultTBB,
                                ErmakovASparMatMultSTL, ErmakovASparMatMultALL>(PPC_SETTINGS_ermakov_a_spar_mat_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ErmakovARunPerfTestSparMatMult::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ErmakovARunPerfTestSparMatMult, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace ermakov_a_spar_mat_mult
