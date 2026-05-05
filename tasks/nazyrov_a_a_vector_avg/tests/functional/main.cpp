#include <gtest/gtest.h>
#include <mpi.h>

#include <numeric>
#include <vector>

#include "../../common/include/common.hpp"
#include "../../mpi/include/ops_mpi.hpp"
#include "../../seq/include/ops_seq.hpp"

namespace nazyrov_a_a_vector_avg {

class VectorAvgTest : public ::testing::Test {
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

TEST_F(VectorAvgTest, SeqTest) {
  InType input = {1, 2, 3, 4, 5};
  OutType expected = 3.0;

  auto task = std::make_shared<VectorAvgSEQ>(input);
  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  EXPECT_DOUBLE_EQ(task->GetOutput(), expected);
}

TEST_F(VectorAvgTest, MpiTest) {
  InType input = {1, 2, 3, 4, 5};
  OutType expected = 3.0;

  auto task = std::make_shared<VectorAvgMPI>(input);
  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    EXPECT_DOUBLE_EQ(task->GetOutput(), expected);
  }
}

}  // namespace nazyrov_a_a_vector_avg
