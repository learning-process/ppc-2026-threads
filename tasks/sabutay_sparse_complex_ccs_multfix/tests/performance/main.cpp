#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "sabutay_sparse_complex_ccs_multfix/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_multfix/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sabutay_sparse_complex_ccs_multfix {
namespace {

CCS BuildRandomCcs(int rows, int cols, int seed, int max_per_col) {
  std::mt19937 gen(static_cast<std::uint32_t>(seed));
  std::uniform_real_distribution<double> re(-3.0, 3.0);
  std::uniform_int_distribution<int> per_col(0, max_per_col);

  CCS m;
  m.row_count = rows;
  m.col_count = cols;
  m.col_start.assign(static_cast<std::size_t>(cols) + 1U, 0);
  m.row_index.clear();
  m.nz.clear();

  for (int jcol = 0; jcol < cols; ++jcol) {
    const int take = per_col(gen);
    std::set<int> pick;
    while (std::cmp_less(pick.size(), static_cast<std::size_t>(take))) {
      const int r = static_cast<int>(gen() % static_cast<std::uint32_t>(rows > 0 ? rows : 1));
      pick.insert(r);
    }
    for (int r : pick) {
      m.row_index.push_back(r);
      m.nz.emplace_back(re(gen), re(gen));
    }
    m.col_start[static_cast<std::size_t>(jcol) + 1U] = static_cast<int>(m.nz.size());
  }
  return m;
}

}  // namespace

class SabutayRunPerfTestThreadsFIX : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    in_ = std::make_tuple(BuildRandomCcs(80, 90, 2027, 7), BuildRandomCcs(90, 70, 4044, 6));
  }

  bool CheckTestOutputData(OutType &out) final {
    const CCS &a = std::get<0>(in_);
    const CCS &b = std::get<1>(in_);
    if (out.row_count != a.row_count || out.col_count != b.col_count) {
      return false;
    }
    if (out.col_count > 0) {
      const int tail = out.col_start[static_cast<std::size_t>(out.col_count)];
      if (out.row_index.size() != static_cast<std::size_t>(tail) || out.nz.size() != static_cast<std::size_t>(tail)) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return in_;
  }

 private:
  InType in_;
};

namespace {

TEST_P(SabutayRunPerfTestThreadsFIX, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SabutaySparseComplexCcsMultFixSEQ>(
    PPC_SETTINGS_sabutay_sparse_complex_ccs_multfix);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SabutayRunPerfTestThreadsFIX::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SabutayRunPerfTestThreadsFIX, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_multfix
