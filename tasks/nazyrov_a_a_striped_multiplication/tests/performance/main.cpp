#include <gtest/gtest.h>
#include <mpi.h>
#include <chrono>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "../../common/include/common.hpp"
#include "../../mpi/include/ops_mpi.hpp"
#include "../../seq/include/ops_seq.hpp"

namespace nazyrov_a_a_striped_multiplication {

std::vector<double> generateMatrix(int n) {
  std::vector<double> mat(n * n);
  for (int i = 0; i < n * n; ++i) {
    mat[i] = rand() % 10;
  }
  return mat;
}

class StripedMultiplicationPerfTest : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
      int argc = 0;
      char** argv = nullptr;
      MPI_Init(&argc, &argv);
    }
  }

  static void TearDownTestSuite() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) {
      MPI_Finalize();
    }
  }
};

TEST_F(StripedMultiplicationPerfTest, SeqPerformance) {
  int n = 100;
  auto A = generateMatrix(n);
  auto B = generateMatrix(n);
  InType input = {A, B};
  
  auto task = std::make_shared<StripedMultiplicationSEQ>(input);
  
  auto start = std::chrono::high_resolution_clock::now();
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
  auto end = std::chrono::high_resolution_clock::now();
  
  double elapsed = std::chrono::duration<double>(end - start).count();
  std::cout << "SEQ time for " << n << "x" << n << ": " << elapsed << "s\n";
  EXPECT_GT(elapsed, 0);
}

TEST_F(StripedMultiplicationPerfTest, MpiPerformance) {
  int n = 100;
  auto A = generateMatrix(n);
  auto B = generateMatrix(n);
  InType input = {A, B};
  
  auto task = std::make_shared<StripedMultiplicationMPI>(input);
  
  auto start = std::chrono::high_resolution_clock::now();
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
  auto end = std::chrono::high_resolution_clock::now();
  
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    double elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "MPI time for " << n << "x" << n << ": " << elapsed << "s\n";
    EXPECT_GT(elapsed, 0);
  }
}

}  // namespace nazyrov_a_a_striped_multiplication
