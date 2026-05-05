#include <gtest/gtest.h>
#include <mpi.h>

#include <chrono>
#include <numeric>
#include <vector>

#include "../../common/include/common.hpp"
#include "../../mpi/include/ops_mpi.hpp"
#include "../../seq/include/ops_seq.hpp"

namespace nazyrov_a_a_vector_avg {

class VectorAvgPerfTest : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
      int argc = 0;
      char **argv = nullptr;
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

TEST_F(VectorAvgPerfTest, SeqPerformance) {
  const int size = 1000000;
  InType input(size, 1);

  auto task = std::make_shared<VectorAvgSEQ>(input);
  auto start = std::chrono::high_resolution_clock::now();
  task->Validation();
  task->PreProcessing();
  task->Run();
  task->PostProcessing();
  auto end = std::chrono::high_resolution_clock::now();

  double elapsed = std::chrono::duration<double>(end - start).count();
  std::cout << "SEQ time: " << elapsed << "s\n";
  EXPECT_GT(elapsed, 0);
}

TEST_F(VectorAvgPerfTest, MpiPerformance) {
  const int size = 1000000;
  InType input(size, 1);

  auto task = std::make_shared<VectorAvgMPI>(input);
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
    std::cout << "MPI time: " << elapsed << "s\n";
    EXPECT_GT(elapsed, 0);
  }
}

}  // namespace nazyrov_a_a_vector_avg
