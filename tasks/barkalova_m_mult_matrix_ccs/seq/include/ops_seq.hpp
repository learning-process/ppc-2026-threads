/*#pragma once

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace barkalova_m_mult_matrix_ccs {

class BarkalobaMMultMatrixCcsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BarkalobaMMultMatrixCcsSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace barkalova_m_mult_matrix_ccs
*/

#pragma once

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace barkalova_m_mult_matrix_ccs {

class BarkalovaMMultMatrixCcsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BarkalovaMMultMatrixCcsSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  // Вспомогательные функции для умножения
  static void TransposeMatrix(const CCSMatrix &src, CCSMatrix &dst);
  static void MultiplyMatrices(const CCSMatrix &a, const CCSMatrix &b, CCSMatrix &c);

  // Порог для отбрасывания близких к нулю значений
  static constexpr double kEpsilon = 1e-12;
};

}  // namespace barkalova_m_mult_matrix_ccs
