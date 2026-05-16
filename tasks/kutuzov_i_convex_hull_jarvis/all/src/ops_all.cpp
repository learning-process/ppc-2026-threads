#include "kutuzov_i_convex_hull_jarvis/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"

namespace kutuzov_i_convex_hull_jarvis {

double KutuzovITestConvexHullALL::s_curr_x = 0.0;
double KutuzovITestConvexHullALL::s_curr_y = 0.0;
double KutuzovITestConvexHullALL::s_epsilon = 1e-9;

double KutuzovITestConvexHullALL::DistanceSquared(double a_x, double a_y, double b_x, double b_y) {
  return ((a_x - b_x) * (a_x - b_x)) + ((a_y - b_y) * (a_y - b_y));
}

double KutuzovITestConvexHullALL::CrossProduct(double o_x, double o_y, double a_x, double a_y, double b_x, double b_y) {
  return ((a_x - o_x) * (b_y - o_y)) - ((a_y - o_y) * (b_x - o_x));
}

bool KutuzovITestConvexHullALL::IsBetterPoint(double cross, double epsilon, double current_x, double current_y,
                                              double i_x, double i_y, double next_x, double next_y) {
  if (cross < -epsilon) {
    return true;
  }
  if (std::abs(cross) < epsilon) {
    return DistanceSquared(current_x, current_y, i_x, i_y) > DistanceSquared(current_x, current_y, next_x, next_y);
  }
  return false;
}

void KutuzovITestConvexHullALL::leftmost_reduce(void *invec, void *inoutvec, int *len, MPI_Datatype *) {
  double *in = static_cast<double *>(invec);
  double *inout = static_cast<double *>(inoutvec);
  for (int i = 0; i < *len; ++i) {
    double x_in = in[3 * i], y_in = in[3 * i + 1];
    double x_io = inout[3 * i], y_io = inout[3 * i + 1];
    if (x_in < x_io || (x_in == x_io && y_in < y_io)) {
      inout[3 * i] = x_in;
      inout[3 * i + 1] = y_in;
      inout[3 * i + 2] = in[3 * i + 2];
    }
  }
}

void KutuzovITestConvexHullALL::next_reduce(void *invec, void *inoutvec, int *len, MPI_Datatype *) {
  double *in = static_cast<double *>(invec);
  double *inout = static_cast<double *>(inoutvec);
  for (int i = 0; i < *len; ++i) {
    double ax = inout[3 * i], ay = inout[3 * i + 1];
    double bx = in[3 * i], by = in[3 * i + 1];
    double cross = CrossProduct(s_curr_x, s_curr_y, ax, ay, bx, by);
    if (IsBetterPoint(cross, s_epsilon, s_curr_x, s_curr_y, bx, by, ax, ay)) {
      inout[3 * i] = bx;
      inout[3 * i + 1] = by;
      inout[3 * i + 2] = in[3 * i + 2];
    }
  }
}

void KutuzovITestConvexHullALL::build_types() {
  MPI_Type_contiguous(3, MPI_DOUBLE, &type_leftmost_);
  MPI_Type_commit(&type_leftmost_);
  type_next_ = type_leftmost_;
}

KutuzovITestConvexHullALL::KutuzovITestConvexHullALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool KutuzovITestConvexHullALL::ValidationImpl() {
  return true;
}

bool KutuzovITestConvexHullALL::PreProcessingImpl() {
  int flag;
  MPI_Initialized(&flag);
  if (!flag) {
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, nullptr);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);

  build_types();

  MPI_Op_create((MPI_User_function *)leftmost_reduce, 1, &op_leftmost_);
  MPI_Op_create((MPI_User_function *)next_reduce, 1, &op_next_);

  auto &input = GetInput();
  int N = static_cast<int>(input.size());

  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank_ != 0) {
    input.resize(N);
  }

  std::vector<double> coords(2 * N);
  if (rank_ == 0) {
    for (int i = 0; i < N; ++i) {
      coords[2 * i] = std::get<0>(input[i]);
      coords[2 * i + 1] = std::get<1>(input[i]);
    }
  }
  MPI_Bcast(coords.data(), 2 * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank_ != 0) {
    for (int i = 0; i < N; ++i) {
      input[i] = {coords[2 * i], coords[2 * i + 1]};
    }
  }

  return true;
}

bool KutuzovITestConvexHullALL::RunImpl() {
  auto &input = GetInput();
  const size_t N = input.size();

  if (N < 3) {
    GetOutput() = input;
    return true;
  }

  size_t start = (N * rank_) / size_;
  size_t end = (N * (rank_ + 1)) / size_;

  double local_lm[3], global_lm[3];
  local_lm[0] = std::get<0>(input[0]);
  local_lm[1] = std::get<1>(input[0]);
  local_lm[2] = 0.0;

#pragma omp parallel default(none) shared(input, start, end, local_lm)
  {
    int li = 0;
    double lx = std::get<0>(input[li]);
    double ly = std::get<1>(input[li]);

#pragma omp for nowait
    for (size_t i = start; i < end; ++i) {
      double x = std::get<0>(input[i]);
      double y = std::get<1>(input[i]);
      if (x < lx || (x == lx && y < ly)) {
        li = i;
        lx = x;
        ly = y;
      }
    }
#pragma omp critical
    {
      if (lx < local_lm[0] || (lx == local_lm[0] && ly < local_lm[1])) {
        local_lm[0] = lx;
        local_lm[1] = ly;
        local_lm[2] = static_cast<double>(li);
      }
    }
  }

  MPI_Allreduce(local_lm, global_lm, 1, type_leftmost_, op_leftmost_, MPI_COMM_WORLD);
  size_t leftmost = static_cast<size_t>(global_lm[2]);

  size_t current = leftmost;
  double current_x = std::get<0>(input[current]);
  double current_y = std::get<1>(input[current]);
  const double epsilon = 1e-9;

  auto &output = GetOutput();
  output.clear();

  while (true) {
    output.push_back(input[current]);

    size_t next = (current + 1) % N;
    double next_x = std::get<0>(input[next]);
    double next_y = std::get<1>(input[next]);

#pragma omp parallel default(none) \
    shared(input, N, current, current_x, current_y, start, end, epsilon, next, next_x, next_y)
    {
      size_t local_next = (current + 1) % N;
      double local_next_x = std::get<0>(input[local_next]);
      double local_next_y = std::get<1>(input[local_next]);

#pragma omp for nowait
      for (size_t i = start; i < end; ++i) {
        if (i == current) {
          continue;
        }
        double i_x = std::get<0>(input[i]);
        double i_y = std::get<1>(input[i]);
        double cross = CrossProduct(current_x, current_y, local_next_x, local_next_y, i_x, i_y);
        if (IsBetterPoint(cross, epsilon, current_x, current_y, i_x, i_y, local_next_x, local_next_y)) {
          local_next = i;
          local_next_x = i_x;
          local_next_y = i_y;
        }
      }

#pragma omp critical
      {
        double cross = CrossProduct(current_x, current_y, next_x, next_y, local_next_x, local_next_y);
        if (IsBetterPoint(cross, epsilon, current_x, current_y, local_next_x, local_next_y, next_x, next_y)) {
          next = local_next;
          next_x = local_next_x;
          next_y = local_next_y;
        }
      }
    }

    s_curr_x = current_x;
    s_curr_y = current_y;
    s_epsilon = epsilon;

    double local_cand[3] = {next_x, next_y, static_cast<double>(next)};
    double global_cand[3];

    MPI_Allreduce(local_cand, global_cand, 1, type_next_, op_next_, MPI_COMM_WORLD);

    current = static_cast<size_t>(global_cand[2]);
    current_x = global_cand[0];
    current_y = global_cand[1];

    if (current == leftmost) {
      break;
    }
  }

  return true;
}

bool KutuzovITestConvexHullALL::PostProcessingImpl() {
  if (op_leftmost_ != MPI_OP_NULL) {
    MPI_Op_free(&op_leftmost_);
  }
  if (op_next_ != MPI_OP_NULL) {
    MPI_Op_free(&op_next_);
  }
  if (type_leftmost_ != MPI_DATATYPE_NULL) {
    MPI_Type_free(&type_leftmost_);
  }
  return true;
}

}  // namespace kutuzov_i_convex_hull_jarvis
