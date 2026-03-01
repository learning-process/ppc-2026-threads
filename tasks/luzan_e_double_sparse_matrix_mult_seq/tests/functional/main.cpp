#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "luzan_e_double_sparse_matrix_mult_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace luzan_e_double_sparse_matrix_mult_seq {

class LuzanEDoubleSparseMatrixMultSeqestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "12" + std::to_string(std::get<0>(test_param)[0]);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string file_name = std::get<0>(params);
    // std::cout << file_name << std::endl;
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_luzan_e_double_sparse_matrix_mult_seq), file_name);
    std::ifstream test_file(abs_path);
    if (!test_file) {
      throw std::runtime_error("Cannot open task file");
    }

    Sparse_matrix A = getFromFile(test_file);
    Sparse_matrix B = getFromFile(test_file);
    test_file.close();

    input_data_ = std::make_tuple(A, B);

    file_name = "ans_" + std::get<0>(params);
    abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_luzan_e_double_sparse_matrix_mult_seq), file_name);
    std::ifstream ans_file(abs_path);
    if (!ans_file) {
      throw std::runtime_error("Cannot open asn file");
    }
    output_data_.getsparsedMatrixFromFile(ans_file);
    ans_file.close();
  }

  bool CheckTestOutputData(OutType &output_data) final {
    std::cout << "\n\noutSize: " << output_data.getVal().size() << "TRUE: " << output_data_.getVal().size() << "\n\n";
    return (output_data == output_data_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType output_data_;
};

namespace {

TEST_P(LuzanEDoubleSparseMatrixMultSeqestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 1> kTestParam = {std::make_tuple("test_1.txt")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<LuzanEDoubleSparseMatrixMultSeq, InType>(
    kTestParam, PPC_SETTINGS_luzan_e_double_sparse_matrix_mult_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    LuzanEDoubleSparseMatrixMultSeqestsThreads::PrintFuncTestName<LuzanEDoubleSparseMatrixMultSeqestsThreads>;

INSTANTIATE_TEST_SUITE_P(FuncTests, LuzanEDoubleSparseMatrixMultSeqestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace luzan_e_double_sparse_matrix_mult_seq
