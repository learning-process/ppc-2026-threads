#include <gtest/gtest.h>
#include <mpi.h>
#include <vector>
#include <cstdlib>
#include <cmath>

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

bool matricesEqual(const std::vector<double>& C1, const std::vector<double>& C2, double eps = 1e-9) {
  if (C1.size() != C2.size()) return false;
  for (size_t i = 0; i < C1.size(); ++i) {
    if (std::abs(C1[i] - C2[i]) > eps) return false;
  }
  return true;
}

class StripedMultiplicationTest : public ::testing::Test {
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

TEST_F(StripedMultiplicationTest, SeqTestSmall) {
  int n = 3;
  auto A = generateMatrix(n);
  auto B = generateMatrix(n);
  InType input = {A, B};
  
  auto task = std::make_shared<StripedMultiplicationSEQ>(input);
  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  EXPECT_EQ(task->GetOutput().size(), n * n);
}

TEST_F(StripedMultiplicationTest, MpiTestSmall) {
  int n = 3;
  auto A = generateMatrix(n);
  auto B = generateMatrix(n);
  InType input = {A, B};
  
  auto task = std::make_shared<StripedMultiplicationMPI>(input);
  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    EXPECT_EQ(task->GetOutput().size(), n * n);
  }
}

TEST_F(StripedMultiplicationTest, SeqVsMpi) {
  int n = 4;
  auto A = generateMatrix(n);
  auto B = generateMatrix(n);
  InType input = {A, B};
  
  auto task_seq = std::make_shared<StripedMultiplicationSEQ>(input);
  task_seq->Validation();
  task_seq->PreProcessing();
  task_seq->Run();
  task_seq->PostProcessing();
  auto C_seq = task_seq->GetOutput();
  
  auto task_mpi = std::make_shared<StripedMultiplicationMPI>(input);
  task_mpi->Validation();
  task_mpi->PreProcessing();
  task_mpi->Run();
  task_mpi->PostProcessing();
  
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    auto C_mpi = task_mpi->GetOutput();
    EXPECT_TRUE(matricesEqual(C_seq, C_mpi));
  }
}

}  // namespace nazyrov_a_a_striped_multiplication
