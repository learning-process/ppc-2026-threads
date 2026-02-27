#include "kichanova_k_lin_system_by_conjug_grad/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>
#include <iostream>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

KichanovaKLinSystemByConjugGradSEQ::KichanovaKLinSystemByConjugGradSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType()  ;
}

bool KichanovaKLinSystemByConjugGradSEQ::ValidationImpl() {
  const InType& input_data = GetInput();
  if (input_data.n <= 0) return false;
  if (input_data.A.size() != static_cast<size_t>(input_data.n * input_data.n)) return false;
  if (input_data.b.size() != static_cast<size_t>(input_data.n)) return false;
  
  return true;
}

bool KichanovaKLinSystemByConjugGradSEQ::PreProcessingImpl() {
  GetOutput().assign(GetInput().n, 0.0);
  return true;
}

bool KichanovaKLinSystemByConjugGradSEQ::RunImpl() {
  const InType& input_data = GetInput();
  OutType& x = GetOutput();
  
  int n = input_data.n;
  if (n == 0) return false;
  
  const std::vector<double>& A = input_data.A;
  const std::vector<double>& b = input_data.b;
  double epsilon = input_data.epsilon;

  std::vector<double> r(n);
  std::vector<double> p(n);
  std::vector<double> Ap(n);
  
  double rr_old = 0.0;
  for (int i = 0; i < n; i++) {
    r[i] = b[i];
    p[i] = r[i];
    rr_old += r[i] * r[i];
  }
  
  double residual_norm = std::sqrt(rr_old);
  if (residual_norm < epsilon) {
    return true;
  }
  
  int max_iter = n * 1000;
  for (int iter = 0; iter < max_iter; iter++) {
    for (int i = 0; i < n; i++) {
      double sum = 0.0;
      const double* A_row = &A[i * n];
      for (int j = 0; j < n; j++) {
        sum += A_row[j] * p[j];
      }
      Ap[i] = sum;
    }
    
    double pAp = 0.0;
    for (int i = 0; i < n; i++) {
      pAp += p[i] * Ap[i];
    }
    
    if (std::abs(pAp) < 1e-30) break;
    
    double alpha = rr_old / pAp;
    
    for (int i = 0; i < n; i++) {
      x[i] += alpha * p[i];
    }
    
    double rr_new = 0.0;
    for (int i = 0; i < n; i++) {
      r[i] -= alpha * Ap[i];
      rr_new += r[i] * r[i];
    }
    
    residual_norm = std::sqrt(rr_new);
    if (residual_norm < epsilon) break;
    
    double beta = rr_new / rr_old;
    
    for (int i = 0; i < n; i++) {
      p[i] = r[i] + beta * p[i];
    }
    
    rr_old = rr_new;
  }
  
  return true;
}

bool KichanovaKLinSystemByConjugGradSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kichanova_k_lin_system_by_conjug_grad
