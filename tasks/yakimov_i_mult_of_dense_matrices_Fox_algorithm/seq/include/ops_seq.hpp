#pragma once

#include <string>

#include "task/include/task.hpp"
#include "yakimov_i_mult_of_dense_matrices_Fox_algorithm/common/include/common.hpp"

namespace yakimov_i_mult_of_dense_matrices_Fox_algorithm {

class YakimovIMultOfDenseMatricesFoxAlgorithmSEQ : public BaseTask {
public:
    static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
        return ppc::task::TypeOfTask::kSEQ;
    }
    explicit YakimovIMultOfDenseMatricesFoxAlgorithmSEQ(const InType& in);

private:
    bool ValidationImpl() override;
    bool PreProcessingImpl() override;
    bool RunImpl() override;
    bool PostProcessingImpl() override;

    DenseMatrix matrix_a_;
    DenseMatrix matrix_b_;
    DenseMatrix result_matrix_;
    std::string matrix_a_filename_;
    std::string matrix_b_filename_;
    int block_size_ = 0;
};

}  // namespace yakimov_i_mult_of_dense_matrices_Fox_algorithm