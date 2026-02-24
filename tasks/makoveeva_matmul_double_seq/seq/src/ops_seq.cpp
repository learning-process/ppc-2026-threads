#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

namespace makoveeva_matmul_double_seq {

MatmulDoubleSeqTask::MatmulDoubleSeqTask(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  
  n_ = std::get<0>(in);
  A_ = std::get<1>(in);
  B_ = std::get<2>(in);
  
  C_.assign(n_ * n_, 0.0);
  GetOutput() = C_;  
}

bool MatmulDoubleSeqTask::ValidationImpl() {

  if (n_ <= 0) return false;
  
  if (A_.size() != n_ * n_ || B_.size() != n_ * n_) return false;
  
  return true;
}

bool MatmulDoubleSeqTask::PreProcessingImpl() {
  return true;
}

bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) return false;
  
 
  int block_size = std::max(1, static_cast<int>(std::sqrt(n_)));
  int num_blocks = (n_ + block_size - 1) / block_size;
  
 
  C_.assign(n_ * n_, 0.0);
  
  for (int ib = 0; ib < num_blocks; ib++) {
    for (int jb = 0; jb < num_blocks; jb++) {
      for (int kb = 0; kb < num_blocks; kb++) {
        int i_start = ib * block_size;
        int i_end = std::min(i_start + block_size, static_cast<int>(n_));
        int j_start = jb * block_size;
        int j_end = std::min(j_start + block_size, static_cast<int>(n_));
        int k_start = kb * block_size;
        int k_end = std::min(k_start + block_size, static_cast<int>(n_));
        
        for (int i = i_start; i < i_end; i++) {
          for (int j = j_start; j < j_end; j++) {
            double sum = 0.0;
            for (int k = k_start; k < k_end; k++) {
              sum += A_[i * n_ + k] * B_[k * n_ + j];
            }
            C_[i * n_ + j] += sum;
          }
        }
      }
    }
  }
  
  GetOutput() = C_;
  
  return true;
}

bool MatmulDoubleSeqTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_seq