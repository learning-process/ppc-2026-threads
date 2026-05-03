#pragma once

#include "kotelnikova_a_double_matr_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kotelnikova_a_double_matr_mult {

class KotelnikovaATaskOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KotelnikovaATaskOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static SparseMatrixCCS MultiplyMatrices(const SparseMatrixCCS &a, const SparseMatrixCCS &b);
  static bool IsMatrixValid(const SparseMatrixCCS &matrix);
};

}  // namespace kotelnikova_a_double_matr_mult
