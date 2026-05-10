#include "../include/ops_seq.hpp"

namespace nazyrov_a_a_striped_multiplication {

StripedMultiplicationSEQ::StripedMultiplicationSEQ(const InType &in) : BaseTask() {
  GetInput() = in;
}

bool StripedMultiplicationSEQ::ValidationImpl() {
  return true;
}

bool StripedMultiplicationSEQ::PreProcessingImpl() {
  auto input = GetInput();
  A_ = input.first;
  B_ = input.second;
  n_ = static_cast<int>(std::sqrt(A_.size()));
  C_.resize(n_ * n_, 0.0);
  return true;
}

bool StripedMultiplicationSEQ::RunImpl() {
  for (int i = 0; i < n_; ++i) {
    for (int k = 0; k < n_; ++k) {
      double aik = A_[i * n_ + k];
      if (aik != 0.0) {
        for (int j = 0; j < n_; ++j) {
          C_[i * n_ + j] += aik * B_[k * n_ + j];
        }
      }
    }
  }
  return true;
}

bool StripedMultiplicationSEQ::PostProcessingImpl() {
  GetOutput() = C_;
  return true;
}

}  // namespace nazyrov_a_a_striped_multiplication