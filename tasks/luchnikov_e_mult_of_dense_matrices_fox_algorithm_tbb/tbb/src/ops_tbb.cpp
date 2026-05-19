#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_tbb/tbb/include/ops_tbb.hpp"
#include <algorithm>
#include <cstddef>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_tbb {

namespace {

void FillMatrixWithOnes(DenseMatrix& mat, int n) {
    mat.rows = n;
    mat.cols = n;
    mat.values.assign(static_cast<std::size_t>(n) * n, 1.0);
}

void FillMatrixWithZeros(DenseMatrix& mat, int r, int c) {
    mat.rows = r;
    mat.cols = c;
    mat.values.assign(static_cast<std::size_t>(r) * c, 0.0);
}

void MultiplyBlockSequential(const DenseMatrix& a, const DenseMatrix& b, DenseMatrix& res,
                             int row_off, int col_off, int blk, int a_col_shift, int b_row_shift) {
    for (int i = 0; i < blk; ++i) {
        for (int j = 0; j < blk; ++j) {
            double acc = 0.0;
            for (int k = 0; k < blk; ++k) {
                acc += a.At(row_off + i, a_col_shift + k) * b.At(b_row_shift + k, col_off + j);
            }
            res.At(row_off + i, col_off + j) += acc;
        }
    }
}

struct FoxStageBody {
    const DenseMatrix* a{};
    const DenseMatrix* b{};
    DenseMatrix* res{};
    int blk{};
    int stage{};
    int num_stages{};

    void operator()(const tbb::blocked_range<int>& range) const {
        for (int i = range.begin(); i < range.end(); ++i) {
            for (int j = 0; j < num_stages; ++j) {
                int broadcast_idx = (i + stage) % num_stages;
                MultiplyBlockSequential(*a, *b, *res,
                                        i * blk, j * blk, blk,
                                        broadcast_idx * blk, broadcast_idx * blk);
            }
        }
    }
};

int DetermineBlockSize(int n) {
    int blk = n / 4;
    if (blk <= 0) { blk = 1;
}
    return std::min(blk, 128);
}

void ExecuteFoxAlgorithmTBB(const DenseMatrix& a, const DenseMatrix& b, DenseMatrix& res, int blk) {
    if (!a.IsSquare() || !b.IsSquare() || a.rows != b.rows || blk <= 0 || a.rows % blk != 0)