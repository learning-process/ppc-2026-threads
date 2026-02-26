#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zorin_d_strassen_alg_matrix_seq/common/include/common.hpp"
#include "zorin_d_strassen_alg_matrix_seq/seq/include/ops_seq.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

std::vector<double> MakeOnes(std::size_t n) {
  return std::vector<double>(n * n, 1.0);
}

class ZorinDRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    constexpr std::size_t n = 256;
    input_.n = n;
    input_.A = MakeOnes(n);
    input_.B = MakeOnes(n);
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

TEST_P(ZorinDRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZorinDStrassenAlgMatrixSEQ>(PPC_SETTINGS_zorin_d_strassen_alg_matrix_seq);

const auto kValues = ppc::util::TupleToGTestValues(kPerfTasks);
const auto kName = ZorinDRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(StrassenPerf, ZorinDRunPerfTests, kValues, kName);

}  // namespace

}  // namespace zorin_d_strassen_alg_matrix_seq
