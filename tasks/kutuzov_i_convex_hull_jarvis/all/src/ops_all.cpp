#include "kutuzov_i_convex_hull_jarvis/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <array>
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

void KutuzovITestConvexHullALL::LeftmostReduce(void *invec, void *inoutvec, const int *len, MPI_Datatype * /*unused*/) {
  auto *in = static_cast<double *>(invec);
  auto *inout = static_cast<double *>(inoutvec);
  for (int i = 0; i < *len; ++i) {
    const auto idx = static_cast<ptrdiff_t>(3) * static_cast<ptrdiff_t>(i);
    double x_in = in[idx];
    double y_in = in[idx + 1];
    double x_io = inout[idx];
    double y_io = inout[idx + 1];
    if (x_in < x_io || (x_in == x_io && y_in < y_io)) {
      inout[idx] = x_in;
      inout[idx + 1] = y_in;
      inout[idx + 2] = in[idx + 2];
    }
  }
}

void KutuzovITestConvexHullALL::NextReduce(void *invec, void *inoutvec, const int *len, MPI_Datatype * /*unused*/) {
  auto *in = static_cast<double *>(invec);
  auto *inout = static_cast<double *>(inoutvec);
  for (int i = 0; i < *len; ++i) {
    const auto idx = static_cast<ptrdiff_t>(3) * static_cast<ptrdiff_t>(i);
    double ax = inout[idx];
    double ay = inout[idx + 1];
    double bx = in[idx];
    double by = in[idx + 1];
    double cross = CrossProduct(s_curr_x, s_curr_y, ax, ay, bx, by); 
    if (IsBetterPoint(cross, s_epsilon, s_curr_x, s_curr_y, bx, by, ax, ay)) {
      inout[idx] = bx;
      inout[idx + 1] = by;
      inout[idx + 2] = in[idx + 2];
    }
  }
}

void KutuzovITestConvexHullALL::BuildTypes() {
  MPI_Type_contiguous(3, MPI_DOUBLE, &type_leftmost_);
  MPI_Type_commit(&type_leftmost_);
  type_next_ = type_leftmost_;
}

KutuzovITestConvexHullALL::KutuzovITestConvexHullALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

size_t KutuzovITestConvexHullALL::FindLeftmostPoint(const InType &input, size_t start, size_t end) {
  std::array<double, 3> local_lm{};
  std::array<double, 3> global_lm{};
  local_lm[0] = std::get<0>(input[0]);
  local_lm[1] = std::get<1>(input[0]);
  local_lm[2] = 0.0;

#pragma omp parallel default(none) shared(input, start, end, local_lm)
  {
    size_t li = 0;
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

  MPI_Allreduce(local_lm.data(), global_lm.data(), 1, type_leftmost_, op_leftmost_, MPI_COMM_WORLD);
  return static_cast<size_t>(global_lm[2]);
}

size_t KutuzovITestConvexHullALL::FindNextPoint(const InType &input, size_t n, size_t current, double current_x,
                                                double current_y, double epsilon, size_t start, size_t end) {
  size_t next = (current + 1) % n;
  double next_x = std::get<0>(input[next]);
  double next_y = std::get<1>(input[next]);

#pragma omp parallel default(none) \
    shared(input, n, current, current_x, current_y, start, end, epsilon, next, next_x, next_y)
  {
    size_t local_next = (current + 1) % n;
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

  std::array<double, 3> local_cand{next_x, next_y, static_cast<double>(next)};
  std::array<double, 3> global_cand{};

  MPI_Allreduce(local_cand.data(), global_cand.data(), 1, type_next_, op_next_, MPI_COMM_WORLD);

  return static_cast<size_t>(global_cand[2]);
}

bool KutuzovITestConvexHullALL::ValidationImpl() {
  return true;
}

bool KutuzovITestConvexHullALL::PreProcessingImpl() {
  int flag = 0;
  MPI_Initialized(&flag);
  if (flag == 0) {
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, nullptr);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);

  BuildTypes();

  MPI_Op_create(reinterpret_cast<MPI_User_function *>(LeftmostReduce), 1, &op_leftmost_);
  MPI_Op_create(reinterpret_cast<MPI_User_function *>(NextReduce), 1, &op_next_);

  auto &input = GetInput();
  auto n = static_cast<int>(input.size());

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank_ != 0) {
    input.resize(static_cast<size_t>(n));
  }

  std::vector<double> coords(static_cast<size_t>(2) * static_cast<size_t>(n));
  if (rank_ == 0) {
    for (int i = 0; i < n; ++i) {
      coords[(2 * static_cast<size_t>(i))] = std::get<0>(input[static_cast<size_t>(i)]);
      coords[(2 * static_cast<size_t>(i)) + 1] = std::get<1>(input[static_cast<size_t>(i)]);
    }
  }
  MPI_Bcast(coords.data(), static_cast<int>(coords.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank_ != 0) {
    for (int i = 0; i < n; ++i) {
      const auto idx = static_cast<size_t>(i);
      input[idx] = {coords[2 * idx], coords[(2 * idx) + 1]};
    }
  }

  return true;
}

bool KutuzovITestConvexHullALL::RunImpl() {
  auto &input = GetInput();
  const auto n = static_cast<size_t>(input.size());

  if (n < 3) {
    GetOutput() = input;
    return true;
  }

  const size_t start = (n * static_cast<size_t>(rank_)) / static_cast<size_t>(size_);
  const size_t end = (n * (static_cast<size_t>(rank_) + 1)) / static_cast<size_t>(size_);
  const double epsilon = 1e-9;

  const size_t leftmost = FindLeftmostPoint(input, start, end);

  size_t current = leftmost;
  double current_x = std::get<0>(input[current]);
  double current_y = std::get<1>(input[current]);

  auto &output = GetOutput();
  output.clear();

  while (true) {
    output.push_back(input[current]);

    const size_t next = FindNextPoint(input, n, current, current_x, current_y, epsilon, start, end);

    if (next == leftmost) {
      break;
    }

    current = next;
    current_x = std::get<0>(input[current]);
    current_y = std::get<1>(input[current]);
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
