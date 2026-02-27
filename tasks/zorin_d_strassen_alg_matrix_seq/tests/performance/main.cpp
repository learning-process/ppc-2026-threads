#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zorin_d_strassen_alg_matrix_seq/common/include/common.hpp"
#include "zorin_d_strassen_alg_matrix_seq/seq/include/ops_seq.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

std::vector<double> make_ones(std::size_t n) {
  return std::vector<double>(n * n, 1.0);
}

class ZorinDRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    constexpr std::size_t k_n = 256;
    input_.n = k_n;
    input_.a = make_ones(k_n);
    input_.b = make_ones(k_n);
  }

  InType GetTestInputData() final {
    return input_;
  }

  bool CheckTestOutputData(OutType &out) final {
    return out.size() == input_.n * input_.n;
  }

 private:
  InType input_{};
};

TEST_P(ZorinDRunPerfTests, ZorinDSEQStrassenRunPerfModes) {
  ExecuteTest(GetParam());
}

const auto k_perf_tasks =
    ppc::util::MakeAllPerfTasks<InType, ZorinDStrassenAlgMatrixSEQ>(PPC_SETTINGS_zorin_d_strassen_alg_matrix_seq);

const auto k_values = ppc::util::TupleToGTestValues(k_perf_tasks);
const auto k_name = ZorinDRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(StrassenPerf, ZorinDRunPerfTests, k_values, k_name);

}  // namespace

}  // namespace zorin_d_strassen_alg_matrix_seq
