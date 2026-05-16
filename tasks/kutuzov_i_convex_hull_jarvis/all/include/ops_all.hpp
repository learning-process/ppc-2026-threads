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

  int rank_ = 0;
  int size_ = 1;
  MPI_Comm comm_ = MPI_COMM_WORLD;
  MPI_Op op_leftmost_ = MPI_OP_NULL;
  MPI_Op op_next_ = MPI_OP_NULL;
  MPI_Datatype type_leftmost_ = MPI_DATATYPE_NULL;
  MPI_Datatype type_next_ = MPI_DATATYPE_NULL;

  static double s_curr_x, s_curr_y, s_epsilon;

  void build_types();

  static void leftmost_reduce(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype);
  static void next_reduce(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype);
};

}  // namespace kutuzov_i_convex_hull_jarvis
