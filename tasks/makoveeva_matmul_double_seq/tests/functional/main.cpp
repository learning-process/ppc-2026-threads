#include <gtest/gtest.h>

#include <vector>
#include <cmath>
#include <tuple>

#include "makoveeva_matmul_double_seq/common/include/common.hpp"
#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

namespace makoveeva_matmul_double_seq {

std::vector<double> SimpleMultiply(const std::vector<double>& A, 
                                   const std::vector<double>& B, 
                                   size_t n) {
  std::vector<double> C(n * n, 0.0);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      double sum = 0.0;
      for (size_t k = 0; k < n; k++) {
        sum += A[i * n + k] * B[k * n + j];
      }
      C[i * n + j] = sum;
    }
  }
  return C;
}


TEST(MatmulDoubleTest, Test2x2) {
  size_t n = 2;
  

  std::vector<double> A = {1.0, 2.0, 3.0, 4.0};

  std::vector<double> B = {5.0, 6.0, 7.0, 8.0};
  std::vector<double> expected = {19.0, 22.0, 43.0, 50.0};
  
  auto input = std::make_tuple(n, A, B);
  MatmulDoubleSeqTask task(input);
  
  ASSERT_TRUE(task.ValidationImpl());
  ASSERT_TRUE(task.PreProcessingImpl());
  ASSERT_TRUE(task.RunImpl());
  ASSERT_TRUE(task.PostProcessingImpl());
  
  auto result = task.GetResult();
  
  for (size_t i = 0; i < result.size(); i++) {
    EXPECT_NEAR(result[i], expected[i], 1e-10);
  }
}

TEST(MatmulDoubleTest, TestVariousSizes) {
  std::vector<size_t> sizes = {2, 3, 4, 5, 8};
  
  for (size_t n : sizes) {
    std::vector<double> A(n * n);
    std::vector<double> B(n * n);
    
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        A[i * n + j] = static_cast<double>(i + j + 1);
        B[i * n + j] = static_cast<double>(i * n + j + 1);
      }
    }
    
    auto expected = SimpleMultiply(A, B, n);
    auto input = std::make_tuple(n, A, B);
    
    MatmulDoubleSeqTask task(input);
    
    ASSERT_TRUE(task.ValidationImpl()) << "Failed validation for n=" << n;
    ASSERT_TRUE(task.PreProcessingImpl()) << "Failed preprocessing for n=" << n;
    ASSERT_TRUE(task.RunImpl()) << "Failed run for n=" << n;
    ASSERT_TRUE(task.PostProcessingImpl()) << "Failed postprocessing for n=" << n;
    
    auto result = task.GetResult();
    
    for (size_t i = 0; i < result.size(); i++) {
      EXPECT_NEAR(result[i], expected[i], 1e-10) << "Failed at index " << i << " for n=" << n;
    }
  }
}

}  // namespace makoveeva_matmul_double_seq