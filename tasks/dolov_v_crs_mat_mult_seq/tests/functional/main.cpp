#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>
#include <string>

#include "dolov_v_crs_mat_mult_seq/common/include/common.hpp"
#include "dolov_v_crs_mat_mult_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dolov_v_crs_mat_mult_seq {

namespace {

SparseMatrix GenerateRandomCRS(int rows, int cols, double density, uint32_t seed) {
    SparseMatrix matrix;
    matrix.num_rows = rows;
    matrix.num_cols = cols;
    matrix.row_pointers.assign(rows + 1, 0);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_real_distribution<> val_dis(-100.0, 100.0);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (dis(gen) < density) {
                matrix.values.push_back(val_dis(gen));
                matrix.col_indices.push_back(j);
            }
        }
        matrix.row_pointers[i + 1] = static_cast<int>(matrix.values.size());
    }
    return matrix;
}

std::vector<double> DenseMultiply(const SparseMatrix& A, const SparseMatrix& B) {
    std::vector<double> res(A.num_rows * B.num_cols, 0.0);
    for (int i = 0; i < A.num_rows; ++i) {
        for (int j = A.row_pointers[i]; j < A.row_pointers[i + 1]; ++j) {
            int col_a = A.col_indices[j];
            double val_a = A.values[j];
            for (int k = B.row_pointers[col_a]; k < B.row_pointers[col_a + 1]; ++k) {
                res[i * B.num_cols + B.col_indices[k]] += val_a * B.values[k];
            }
        }
    }
    return res;
}

} // namespace

class DolovVCrsMatMultSeqRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "Test_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    int test_id = std::get<0>(std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));
    
    SparseMatrix A, B;
    if (test_id == 1) {
        A = GenerateRandomCRS(10, 10, 0.1, 42);
        B = GenerateRandomCRS(10, 10, 0.1, 43);
    } else if (test_id == 2) {
        A = GenerateRandomCRS(20, 10, 0.2, 44);
        B = GenerateRandomCRS(10, 15, 0.2, 45);
    } else if (test_id == 3) {
        A = GenerateRandomCRS(5, 5, 0.5, 10);
        B.num_rows = 5; B.num_cols = 5; B.row_pointers.assign(6, 0);
    } else if (test_id == 4) {
        A = GenerateRandomCRS(1, 100, 0.5, 123);
        B = GenerateRandomCRS(100, 1, 0.5, 321);
    } else if (test_id == 5) {
        A = GenerateRandomCRS(100, 100, 0.02, 49);
        B = GenerateRandomCRS(100, 100, 0.02, 50);
    } else if (test_id == 6) {
        A.num_rows = 3; A.num_cols = 3;
        A.row_pointers = {0, 3, 3, 3};
        A.col_indices = {0, 1, 2}; A.values = {1, 1, 1};
        B = GenerateRandomCRS(3, 3, 0.5, 99);
    } else {
        A = GenerateRandomCRS(1, 1, 1.0, 1);
        B = GenerateRandomCRS(1, 1, 1.0, 2);
    }

    expected_dense_ = DenseMultiply(A, B);
    rows_res_ = A.num_rows;
    cols_res_ = B.num_cols;
    input_data_ = {A, B};
  }

  bool CheckTestOutputData(OutType &out) final {
    if (out.num_rows != rows_res_ || out.num_cols != cols_res_) return false;
    
    std::vector<double> actual_dense(out.num_rows * out.num_cols, 0.0);
    for (int i = 0; i < out.num_rows; ++i) {
        for (int j = out.row_pointers[i]; j < out.row_pointers[i + 1]; ++j) {
            actual_dense[i * out.num_cols + out.col_indices[j]] = out.values[j];
        }
    }

    for (size_t i = 0; i < actual_dense.size(); ++i) {
        if (std::abs(actual_dense[i] - expected_dense_[i]) > 1e-9) return false;
    }
    return true;
  }

  InType GetTestInputData() final { return input_data_; }

 private:
  InType input_data_;
  std::vector<double> expected_dense_;
  int rows_res_, cols_res_;
};

namespace {

TEST_P(DolovVCrsMatMultSeqRunFuncTestsThreads, RandomSparseMatrices) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(1, "SmallSparse"), 
    std::make_tuple(2, "Rectangular"),
    std::make_tuple(3, "ZeroMatrix"), 
    std::make_tuple(4, "DotProduct"),
    std::make_tuple(5, "LargeSparse"),
    std::make_tuple(6, "DenseRow"),
    std::make_tuple(7, "SingleElement")
};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<DolovVCrsMatMultSeq, InType>(kTestParam, PPC_SETTINGS_dolov_v_crs_mat_mult_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = DolovVCrsMatMultSeqRunFuncTestsThreads::PrintFuncTestName<DolovVCrsMatMultSeqRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(CRS_Functional_Advanced_Tests, DolovVCrsMatMultSeqRunFuncTestsThreads, kGtestValues, kPerfTestName);

} // namespace
} // namespace dolov_v_crs_mat_mult_seq