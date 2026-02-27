#include "timur_a_cannon/seq/include/ops_seq.hpp"

#include <utility>
#include <vector>

#include "timur_a_cannon/common/include/common.hpp"

namespace timur_a_cannon {


struct InputData {
    int block_size;
    const std::vector<std::vector<double>>& matrix_a;
    const std::vector<std::vector<double>>& matrix_b;
};

InputData GetInputData(const InType& input_tuple) {
    return {
        std::get<0>(input_tuple),
        std::get<1>(input_tuple),
        std::get<2>(input_tuple)
    };
}

bool TimurACannonMatrixMultiplication::ValidationImpl() {
    const auto& input = GetInputData(GetInput());
    const auto& matrix_a = input.matrix_a;
    const auto& matrix_b = input.matrix_b;
    const int block_size = input.block_size;

    if (block_size <= 0) {
        return false;
    }

    if (matrix_a.empty()  matrix_b.empty()) {
        return false;
    }

    if (matrix_a.size() != matrix_a[0].size()  matrix_b.size() != matrix_b[0].size()) {
        return false;
    }

    if (matrix_a.size() != matrix_b.size()) {
        return false;
    }

    if (matrix_a.size() % static_cast<size_t>(block_size) != 0) {
        return false;
    }

    return true;
}

bool TimurACannonMatrixMultiplication::PreProcessingImpl() {
    GetOutput().clear();
    return true;
}

void TimurACannonMatrixMultiplication::MultiplyingBlocks(
    std::vector<std::vector<double>>& result_block, 
    const std::vector<std::vector<double>>& block_a, 
    const std::vector<std::vector<double>>& block_b, 
    int block_dimension)
{

    for (int i = 0; i < block_dimension; ++i) {
        for (int j = 0; j < block_dimension; ++j) {
            for (int k = 0; k < block_dimension; ++k) {
                result_block[i][j] += block_a[i][k] * block_b[k][j];
            }
        }
    }
}

void TimurACannonMatrixMultiplication::ShiftBlocksMatrixALeft(
    std::vector<std::vector<std::vector<std::vector<double>>>>& matrix_blocks_a,
    int num_block_columns)
{
    for (int row_block_idx = 0; row_block_idx < num_block_columns; ++row_block_idx) {
        auto temp_block = std::move(matrix_blocks_a[row_block_idx][0]);
        for (int col_block_idx = 1; col_block_idx < num_block_columns; ++col_block_idx) {
            matrix_blocks_a[row_block_idx][col_block_idx - 1] = std::move(matrix_blocks_a[row_block_idx][col_block_idx]);
        }
        matrix_blocks_a[row_block_idx][num_block_columns - 1] = std::move(temp_block);
    }
}

void TimurACannonMatrixMultiplication::ShiftBlocksMatrixBUp(
    std::vector<std::vector<std::vector<std::vector<double>>>>& matrix_blocks_b,
    int num_block_rows)
{
    for (int col_block_idx = 0; col_block_idx < num_block_rows; ++col_block_idx) {
        auto temp_block = std::move(matrix_blocks_b[0][col_block_idx]);
        for (int row_block_idx = 1; row_block_idx < num_block_rows; ++row_block_idx) {
            matrix_blocks_b[row_block_idx - 1][col_block_idx] = std::move(matrix_blocks_b[row_block_idx][col_block_idx]);
        }
        matrix_blocks_b[num_block_rows - 1][col_block_idx] = std::move(temp_block);
    }
}

void TimurACannonMatrixMultiplication::AlgorithmCannon(
    std::vector<std::vector<std::vector<std::vector<double>>>>& matrix_blocks_a,
    std::vector<std::vector<std::vector<std::vector<double>>>>& matrix_blocks_b,
    std::vector<std::vector<std::vector<std::vector<double>>>>& matrix_blocks_c,
    int block_dimension,
    int num_block_dim) 
{
    for (int step = 0; step < num_block_dim; ++step) {
        for (int i = 0; i < num_block_dim; ++i) {
            for (int j = 0; j < num_block_dim; ++j) {
                MultiplyingBlocks(matrix_blocks_c[i][j], matrix_blocks_a[i][j], matrix_blocks_b[i][j], block_dimension);
            }
        }
        if (step < num_block_dim - 1) {
            ShiftBlocksMatrixALeft(matrix_blocks_a, num_block_dim);
            ShiftBlocksMatrixBUp(matrix_blocks_b, num_block_dim);
        }
    }
}

void TimurACannonMatrixMultiplication::FillingResultingMatrix(
    const std::vector<std::vector<std::vector<std::vector<double>>>>& matrix_blocks_c,
    std::vector<std::vector<double>>& result_matrix,
    int block_dimension,
    int num_block_dim)
{
    for (int i = 0; i < num_block_dim; ++i) { 
        for (int j = 0; j < num_block_dim; ++j) { 
            for (int row = 0; row < block_dimension; ++row) { 
                for (int col = 0; col < block_dimension; ++col) { 
                    result_matrix[(i * block_dimension) + row][(j * block_dimension) + col] =
                        matrix_blocks_c[i][j][row][col];
                }
            }
        }
    }
}

bool TimurACannonMatrixMultiplication::RunImpl() {
    const auto& input = GetInputData(GetInput());
    const int block_size = input.block_size;
    const auto& matrix_a = input.matrix_a;
    const auto& matrix_b = input.matrix_b;
    const int matrix_dimension = static_cast<int>(matrix_a.size());
    const int num_block_dim = matrix_dimension / block_size; 
    std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_a(
        num_block_dim,
        std::vector<std::vector<std::vector<double>>>(
            num_block_dim, std::vector<std::vector<double>>(block_size, std::vector<double>(block_size))));

    std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_b(
        num_block_dim,
        std::vector<std::vector<std::vector<double>>>(
            num_block_dim, std::vector<std::vector<double>>(block_size, std::vector<double>(block_size))));

    std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_c(
        num_block_dim,
        std::vector<std::vector<std::vector<double>>>(
            num_block_dim, std::vector<std::vector<double>>(block_size, std::vector<double>(block_size, 0.0))));
    for (int i = 0; i < num_block_dim; ++i) { 
        for (int j = 0; j < num_block_dim; ++j) { 
            int shift_a_col = (j + i) % num_block_dim;
            int shift_b_row = (i + j) % num_block_dim;
            for (int row = 0; row < block_size; ++row) {
                for (int col = 0; col < block_size; ++col) {
                    matrix_blocks_a[i][j][row][col] = matrix_a[(i * block_size) + row][(shift_a_col * block_size) + col];
                    matrix_blocks_b[i][j][row][col] = matrix_b[(shift_b_row * block_size) + row][(j * block_size) + col];
                }
            }
        }
    }

    AlgorithmCannon(matrix_blocks_a, matrix_blocks_b, matrix_blocks_c, block_size, num_block_dim);
    std::vector<std::vector<double>> result_matrix(matrix_dimension, std::vector<double>(matrix_dimension));
    FillingResultingMatrix(matrix_blocks_c, result_matrix, block_size, num_block_dim);
    GetOutput() = std::move(result_matrix);
    return true;
}

bool TimurACannonMatrixMultiplication::PostProcessingImpl() {

    return true;
}

} // namespace timur_a_cannon