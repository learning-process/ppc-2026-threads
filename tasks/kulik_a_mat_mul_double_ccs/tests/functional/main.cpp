#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <ios>
#include <string>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"
#include "kulik_a_mat_mul_double_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kulik_a_mat_mul_double_ccs {

class KulikARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string filename_a = std::get<0>(params) + ".bin";
    std::string filename_b = std::get<1>(params) + ".bin";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_kulik_a_mat_mul_double_ccs, filename_a);
    std::ifstream filestream(abs_path, std::ios::in | std::ios::binary);
    CCS matrix_test_a;
    filestream.read(reinterpret_cast<char *>(&matrix_test_a.n), sizeof(size_t));
    filestream.read(reinterpret_cast<char *>(&matrix_test_a.m), sizeof(size_t));
    filestream.read(reinterpret_cast<char *>(&matrix_test_a.nz), sizeof(size_t));
    matrix_test_a.col_ind.resize(matrix_test_a.m + 1);
    matrix_test_a.row.resize(matrix_test_a.nz);
    matrix_test_a.value.resize(matrix_test_a.nz);
    filestream.read(reinterpret_cast<char *>(matrix_test_a.col_ind.data()),
                    static_cast<std::streamsize>(matrix_test_a.col_ind.size() * sizeof(size_t)));
    filestream.read(reinterpret_cast<char *>(matrix_test_a.row.data()),
                    static_cast<std::streamsize>(matrix_test_a.row.size() * sizeof(size_t)));
    filestream.read(reinterpret_cast<char *>(matrix_test_a.value.data()),
                    static_cast<std::streamsize>(matrix_test_a.value.size() * sizeof(double)));

    filestream.close();
    abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_kulik_a_mat_mul_double_ccs, filename_b);
    std::ifstream filestream2(abs_path, std::ios::in | std::ios::binary);
    CCS matrix_test_b;
    filestream2.read(reinterpret_cast<char *>(&matrix_test_b.n), sizeof(size_t));
    filestream2.read(reinterpret_cast<char *>(&matrix_test_b.m), sizeof(size_t));
    filestream2.read(reinterpret_cast<char *>(&matrix_test_b.nz), sizeof(size_t));
    matrix_test_b.col_ind.resize(matrix_test_b.m + 1);
    matrix_test_b.row.resize(matrix_test_b.nz);
    matrix_test_b.value.resize(matrix_test_b.nz);
    filestream2.read(reinterpret_cast<char *>(matrix_test_b.col_ind.data()),
                     static_cast<std::streamsize>(matrix_test_b.col_ind.size() * sizeof(size_t)));
    filestream2.read(reinterpret_cast<char *>(matrix_test_b.row.data()),
                     static_cast<std::streamsize>(matrix_test_b.row.size() * sizeof(size_t)));
    filestream2.read(reinterpret_cast<char *>(matrix_test_b.value.data()),
                     static_cast<std::streamsize>(matrix_test_b.value.size() * sizeof(double)));

    filestream2.close();
    input_data_ = std::make_tuple(matrix_test_a, matrix_test_b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    std::vector<double> value{1, 35, 42, 5, 27, 32, 64, 28, 2, 20, 70, 36};
    std::vector<size_t> row{0, 3, 5, 3, 5, 1, 3, 5, 0, 1, 3, 5};
    std::vector<size_t> col_ind{0, 1, 3, 5, 8, 9, 12};
    bool f1 = true;
    bool f2 = true;
    bool f3 = true;
    for (size_t i = 0; i < value.size(); ++i) {
      if (std::abs(output_data.value[i] - value[i]) > 1e-12) {
        f1 = false;
      }
    }
    f2 = (output_data.row == row);
    f3 = (output_data.col_ind == col_ind);
    return (f1 && f2 && f3);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(KulikARunFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 1> kTestParam = {std::make_tuple(std::string("matrix_test"), std::string("matrix_test"))};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KulikAMatMulDoubleCcsSEQ, InType>(kTestParam, PPC_SETTINGS_kulik_a_mat_mul_double_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KulikARunFuncTestsThreads::PrintFuncTestName<KulikARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, KulikARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kulik_a_mat_mul_double_ccs
