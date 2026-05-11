#include "safronov_m_multiplication_matrix_blockscheme_cannon/all/include/ops_all.hpp"
#include <cmath>
#include <mpi.h>
#include "oneapi/tbb/blocked_range2d.h"
#include "oneapi/tbb/parallel_for.h"
#include <algorithm>

namespace safronov_m_multiplication_matrix_blocksscheme_cannon {

SafronovMMultiplicationMatrixBlockSchemeCannonALL::SafronovMMultiplicationMatrixBlockSchemeCannonALL(const InType &in) {
    SetTypeOfTask(GetStaticTypeOfTask());
    GetInput() = in;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonALL::ValidationImpl() {
    const auto &in = GetInput();
    int size_block = std::get<0>(in);
    const auto &matrix_a = std::get<1>(in);
    const auto &matrix_b = std::get<2>(in);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int q = static_cast<int>(std::sqrt(size));
    return (size_block > 0) && (!matrix_a.empty() && !matrix_b.empty()) && (matrix_a.size() == matrix_a[0].size()) &&
         (matrix_b.size() == matrix_b[0].size()) && (matrix_a.size() == matrix_b.size()) &&
         (matrix_a.size() % size_block == 0) && (q * q == size) && (matrix_a.size() % q == 0) &&
        (size_block == static_cast<int>(matrix_a.size()) / q);
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonALL::PreProcessingImpl() {
    return true;
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::ParallelMultiplyBlocks(
    const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C, int size_block) {
    
    tbb::parallel_for(tbb::blocked_range2d<int>(0, size_block, 0, size_block),
        [&](const tbb::blocked_range2d<int>& r) {
            for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
                for (int k = 0; k < size_block; ++k) {
                    double temp = A[i * size_block + k];
                    for (int j = r.cols().begin(); j < r.cols().end(); ++j) {
                        C[i * size_block + j] += temp * B[k * size_block + j];
                    }
                }
            }
        });
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::DistributeData(
    int rank, int size, int q, int size_block, std::vector<double>& local_A, std::vector<double>& local_B) {
    
    if (rank == 0) {
        const auto &matrix_a_full = std::get<1>(GetInput());
        const auto &matrix_b_full = std::get<2>(GetInput());

        for (int p = 0; p < size; ++p) {
            int p_row = p / q;
            int p_col = p % q;
            std::vector<double> send_A(size_block * size_block);
            std::vector<double> send_B(size_block * size_block);

            for (int i = 0; i < size_block; ++i) {
                for (int j = 0; j < size_block; ++j) {
                    int a_col = ((p_col + p_row) % q) * size_block + j;
                    int a_row = p_row * size_block + i;
                    send_A[i * size_block + j] = matrix_a_full[a_row][a_col];

                    int b_row = ((p_row + p_col) % q) * size_block + i;
                    int b_col = p_col * size_block + j;
                    send_B[i * size_block + j] = matrix_b_full[b_row][b_col];
                }
            }

            if (p == 0) {
                local_A = std::move(send_A);
                local_B = std::move(send_B);
            } else {
                MPI_Send(send_A.data(), size_block * size_block, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
                MPI_Send(send_B.data(), size_block * size_block, MPI_DOUBLE, p, 1, MPI_COMM_WORLD);
            }
        }
    } else {
        MPI_Recv(local_A.data(), size_block * size_block, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_B.data(), size_block * size_block, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::CannonAlgorithm(
    int rank, int q, int size_block, std::vector<double>& local_A, std::vector<double>& local_B, std::vector<double>& local_C) {
    
    int row = rank / q;
    int col = rank % q;
    int left = row * q + (col - 1 + q) % q;
    int right = row * q + (col + 1) % q;
    int up = ((row - 1 + q) % q) * q + col;
    int down = ((row + 1) % q) * q + col;

    for (int step = 0; step < q; ++step) {
        ParallelMultiplyBlocks(local_A, local_B, local_C, size_block);

        if (step < q - 1) {
            std::vector<double> next_A(size_block * size_block);
            std::vector<double> next_B(size_block * size_block);

            MPI_Sendrecv(local_A.data(), size_block * size_block, MPI_DOUBLE, left, 10,
                         next_A.data(), size_block * size_block, MPI_DOUBLE, right, 10,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            MPI_Sendrecv(local_B.data(), size_block * size_block, MPI_DOUBLE, up, 11,
                         next_B.data(), size_block * size_block, MPI_DOUBLE, down, 11,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            local_A = std::move(next_A);
            local_B = std::move(next_B);
        }
    }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::CollectAndBroadcast(
    int rank, int size, int q, int size_block, const std::vector<double>& local_C) {
    
    int n = q * size_block;
    std::vector<double> flat_result(n * n);

    if (rank == 0) {
        // Копируем свой блок (rank 0)
        for (int i = 0; i < size_block; ++i)
            for (int j = 0; j < size_block; ++j)
                flat_result[i * n + j] = local_C[i * size_block + j];

        // Собираем блоки от остальных
        std::vector<double> recv_buf(size_block * size_block);
        for (int p = 1; p < size; ++p) {
            MPI_Recv(recv_buf.data(), size_block * size_block, MPI_DOUBLE, p, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int p_row = p / q;
            int p_col = p % q;
            for (int i = 0; i < size_block; ++i) {
                for (int j = 0; j < size_block; ++j) {
                    flat_result[(p_row * size_block + i) * n + (p_col * size_block + j)] = recv_buf[i * size_block + j];
                }
            }
        }
    } else {
        MPI_Send(local_C.data(), size_block * size_block, MPI_DOUBLE, 0, 20, MPI_COMM_WORLD);
    }

    // Рассылаем итоговую матрицу всем процессам
    MPI_Bcast(flat_result.data(), n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Сохраняем в GetOutput на всех процессах
    std::vector<std::vector<double>> final_matrix(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            final_matrix[i][j] = flat_result[i * n + j];
        }
    }
    GetOutput() = std::move(final_matrix);
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonALL::RunImpl() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int q = static_cast<int>(std::sqrt(size));
    int size_block = 0;

    if (rank == 0) 
      size_block = std::get<0>(GetInput());
    MPI_Bcast(&size_block, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<double> local_A(size_block * size_block);
    std::vector<double> local_B(size_block * size_block);
    std::vector<double> local_C(size_block * size_block, 0.0);

    // 1. Рассылка начальных данных
    DistributeData(rank, size, q, size_block, local_A, local_B);

    // 2. Алгоритм Кэннона (основные вычисления)
    CannonAlgorithm(rank, q, size_block, local_A, local_B, local_C);

    // 3. Сборка и финальная рассылка всем процессам
    CollectAndBroadcast(rank, size, q, size_block, local_C);

    return true;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonALL::PostProcessingImpl() {
    return true;
}

} // namespace safronov_m_multiplication_matrix_blocksscheme_cannon