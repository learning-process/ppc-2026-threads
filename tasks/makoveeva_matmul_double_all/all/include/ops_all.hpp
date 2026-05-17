#pragma once

#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_all/common/include/common.hpp"
#include "task/include/task.hpp"

namespace makoveeva_matmul_double_all {

class MatmulDoubleAllTask : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit MatmulDoubleAllTask(const InType& in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] const std::vector<double>& GetResult() const { return c_; }

  using BaseTask::GetOutput;
  using BaseTask::GetInput;

 private:
  void choose_implementation();
  void multiply_seq();
  void multiply_omp();
  void multiply_tbb();
  void multiply_stl();

  size_t n_{0};
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> c_;
  
  enum class TechType { kSeq, kOmp, kTbb, kStl };
  TechType selected_tech_{TechType::kSeq};
};

}  // namespace makoveeva_matmul_double_all