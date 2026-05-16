#include <gtest/gtest.h>

#include <mpi.h>
#include <omp.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "sabutay_sparse_complex_ccs_mult_all/all/include/ops_all.hpp"
#include "sabutay_sparse_complex_ccs_mult_all/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult_all/omp/include/ops_omp.hpp"
#include "sabutay_sparse_complex_ccs_mult_all/seq/include/ops_seq.hpp"
#include "sabutay_sparse_complex_ccs_mult_all/stl/include/ops_stl.hpp"
#include "sabutay_sparse_complex_ccs_mult_all/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace sabutay_sparse_complex_ccs_mult_all {
namespace {

void SetPerfAttributesForTask(ppc::task::TypeOfTask kind, ppc::performance::PerfAttr &perf_attrs) {
  if (kind == ppc::task::TypeOfTask::kMPI || kind == ppc::task::TypeOfTask::kALL) {
    const double t0 = ppc::util::GetTimeMPI();
    perf_attrs.current_timer = std::function<double()>([t0]() { return ppc::util::GetTimeMPI() - t0; });
  } else if (kind == ppc::task::TypeOfTask::kOMP) {
    const double t0 = omp_get_wtime();
    perf_attrs.current_timer = std::function<double()>([t0]() { return omp_get_wtime() - t0; });
  } else if (kind == ppc::task::TypeOfTask::kSEQ || kind == ppc::task::TypeOfTask::kSTL ||
             kind == ppc::task::TypeOfTask::kTBB) {
    const auto t0 = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = std::function<double()>([t0]() {
      const auto now = std::chrono::high_resolution_clock::now();
      return std::chrono::duration<double>(now - t0).count();
    });
  } else {
    throw std::runtime_error("The task type is not supported for performance testing.");
  }
}


CCS BuildRandomCcs(int rows, int cols, int seed, int max_per_col) {
  std::mt19937 gen(static_cast<std::uint32_t>(seed));
  std::uniform_real_distribution<double> re((-3.0), (3.0));
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

class SabutayRunPerfTestThreadsALL : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void ExecuteTest(const ppc::util::PerfTestParam<InType, OutType> &perf_test_param) {
    auto task_getter = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTaskGetter)>(perf_test_param);
    auto test_name = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kNameTest)>(perf_test_param);
    auto mode = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(perf_test_param);

    ASSERT_FALSE(test_name.find("unknown") != std::string::npos);
    if (test_name.find("disabled") != std::string::npos) {
      GTEST_SKIP();
    }

    const auto test_env_scope = ppc::util::test::MakePerTestEnvForCurrentGTest(test_name);

    auto task = task_getter(GetTestInputData());
    ppc::performance::Perf perf(task);
    ppc::performance::PerfAttr perf_attr;
    SetPerfAttributesForTask(task->GetDynamicTypeOfTask(), perf_attr);

    if (mode == ppc::performance::PerfResults::TypeOfRunning::kPipeline) {
      perf.PipelineRun(perf_attr);
    } else if (mode == ppc::performance::PerfResults::TypeOfRunning::kTaskRun) {
      perf.TaskRun(perf_attr);
    } else {
      throw std::runtime_error("The type of performance check for the task was not selected.");
    }

    if (ppc::util::GetMPIRank() == 0) {
      perf.PrintPerfStatistic(test_name);
    }

    OutType output_data = task->GetOutput();
    ASSERT_TRUE(CheckTestOutputData(output_data));
  }

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

TEST_P(SabutayRunPerfTestThreadsALL, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SabutaySparseComplexCcsMultAll, SabutaySparseComplexCcsMultOmpFix,
                                SabutaySparseComplexCcsMultFixSEQ, SabutaySparseComplexCcsMultSTL,
                                SabutaySparseComplexCcsMultFixTBB>(PPC_SETTINGS_sabutay_sparse_complex_ccs_mult_all);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SabutayRunPerfTestThreadsALL::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SabutayRunPerfTestThreadsALL, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult_all
