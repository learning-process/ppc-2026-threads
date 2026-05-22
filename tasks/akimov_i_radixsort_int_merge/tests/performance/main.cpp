#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <random>

#include "akimov_i_radixsort_int_merge/all/include/ops_all.hpp"
#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "akimov_i_radixsort_int_merge/omp/include/ops_omp.hpp"
#include "akimov_i_radixsort_int_merge/seq/include/ops_seq.hpp"
#include "akimov_i_radixsort_int_merge/stl/include/ops_stl.hpp"
#include "akimov_i_radixsort_int_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace akimov_i_radixsort_int_merge {

class AkimovIRadixSortIntMergePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  InType input_data_;
  InType expected_sorted_;

  void SetUp() override {
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(-1000, 1000);
    input_data_.resize(kCount_);
    for (int &val : input_data_) {
      val = dist(gen);
    }
    expected_sorted_ = input_data_;
    std::ranges::sort(expected_sorted_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    if (ppc::util::IsUnderMpirun()) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    if (rank != 0) {
      return true;
    }
    return output_data == expected_sorted_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(AkimovIRadixSortIntMergePerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, AkimovIRadixSortIntMergeSEQ, AkimovIRadixSortIntMergeOMP,
                                AkimovIRadixSortIntMergeTBB, AkimovIRadixSortIntMergeSTL, AkimovIRadixSortIntMergeALL>(
        PPC_SETTINGS_akimov_i_radixsort_int_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = AkimovIRadixSortIntMergePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, AkimovIRadixSortIntMergePerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace akimov_i_radixsort_int_merge
