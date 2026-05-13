#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/common/include/common.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/omp/include/ops_omp.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/seq/include/ops_seq.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/stl/include/ops_stl.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/tbb/include/ops_tbb.hpp"

namespace vinyaikina_e_multidimensional_integrals_simpson_method {
class VinyaikinaESimpsonPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    auto xyz_3d = [](const std::vector<double> &x) { return x[0] * x[1] * x[2]; };
    std::vector<std::pair<double, double>> lims = {{0.0, 0.75}, {0.0, 0.75}, {0.0, 0.75}};
    input_ = InType{0.005, lims, xyz_3d};
    etalon_ = 0.022247314453125;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::fabs(output_data - etalon_) <= 1e-3;
  }

  InType GetTestInputData() final {
    return input_;
  }

 private:
  InType input_;
  OutType etalon_ = 0.0;
};

TEST_P(VinyaikinaESimpsonPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}
namespace {

const auto kPerfTaskName = PPC_SETTINGS_vinyaikina_e_multidimensional_integrals_simpson_method;

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, VinyaikinaEMultidimIntegrSimpsonSEQ, VinyaikinaEMultidimIntegrSimpsonOMP,
                                VinyaikinaEMultidimIntegrSimpsonTBB, VinyaikinaEMultidimIntegrSimpsonSTL>(
        kPerfTaskName);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = VinyaikinaESimpsonPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, VinyaikinaESimpsonPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
