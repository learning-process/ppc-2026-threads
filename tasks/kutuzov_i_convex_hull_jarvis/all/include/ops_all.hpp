#pragma once

#include <mpi.h>

#include <cstddef>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kutuzov_i_convex_hull_jarvis {

class KutuzovITestConvexHullALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KutuzovITestConvexHullALL(const InType &in);

 private:
  static double DistanceSquared(double a_x, double a_y, double b_x, double b_y);
  static double CrossProduct(double o_x, double o_y, double a_x, double a_y, double b_x, double b_y);
  static bool IsBetterPoint(double cross, double epsilon, double current_x, double current_y, double i_x, double i_y,
                            double next_x, double next_y);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  size_t FindLeftmostPoint(const InType &input, size_t start, size_t end);
  size_t FindNextPoint(const InType &input, size_t n, size_t current, double current_x, double current_y,
                       double epsilon, size_t start, size_t end);

  int rank_ = 0;
  int size_ = 1;
  MPI_Op op_leftmost_ = MPI_OP_NULL;
  MPI_Op op_next_ = MPI_OP_NULL;
  MPI_Datatype type_leftmost_ = MPI_DATATYPE_NULL;
  MPI_Datatype type_next_ = MPI_DATATYPE_NULL;

  static double s_curr_x;
  static double s_curr_y;
  static double s_epsilon;

  void BuildTypes();

  static void LeftmostReduce(void *invec, void *inoutvec, const int *len, MPI_Datatype * /*unused*/);
  static void NextReduce(void *invec, void *inoutvec, const int *len, MPI_Datatype * /*unused*/);
};

}  // namespace kutuzov_i_convex_hull_jarvis
