#include "lukin_i_ench_contr_lin_hist/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <vector>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"

namespace lukin_i_ench_contr_lin_hist {

LukinITestTaskALL::LukinITestTaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  }

  GetOutput() = OutType(GetInput().size());
}

bool LukinITestTaskALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    return (static_cast<int>(GetInput().size()) != 0);
  }
  return true;
}

bool LukinITestTaskALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    vec_size_ = static_cast<int>(GetInput().size());
  }
  return true;
}

bool LukinITestTaskALL::RunImpl() {
  int proc_count = -1;
  int proc_rank = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

  MPI_Bcast(&vec_size_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  GetOutput().resize(vec_size_);

  std::vector<int> sendcounts(proc_count);
  std::vector<int> offsets(proc_count);

  const int part = vec_size_ / proc_count;
  const int reminder = vec_size_ % proc_count;

  int offset = 0;
  for (int i = 0; i < proc_count; i++) {
    sendcounts[i] = part + (i < reminder ? 1 : 0);
    offsets[i] = offset;
    offset += sendcounts[i];
  }

  int local_size = sendcounts[proc_rank];
  OutType local_vec(local_size);

  MPI_Scatterv(proc_rank == 0 ? GetInput().data() : nullptr, sendcounts.data(), offsets.data(), MPI_UNSIGNED_CHAR,
               local_vec.data(), local_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  unsigned char local_min = 255;
  unsigned char local_max = 0;

#pragma omp parallel for default(none) shared(local_vec, local_size) reduction(min : local_min) \
    reduction(max : local_max)
  for (int i = 0; i < local_size; i++) {
    local_min = std::min(local_min, local_vec[i]);
    local_max = std::max(local_max, local_vec[i]);
  }

  unsigned char min = 255;
  unsigned char max = 0;

  MPI_Allreduce(&local_min, &min, 1, MPI_UNSIGNED_CHAR, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_max, &max, 1, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD);

  if (max == min) {
    if (proc_rank == 0) {
      GetOutput() = GetInput();
    }
    OutType output = GetOutput();
    MPI_Bcast(output.data(), vec_size_, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    GetOutput() = output;
    return true;
  }

  float scale = 255.0F / static_cast<float>(max - min);

  OutType output(local_size);

#pragma omp parallel for default(none) shared(local_vec, output, min, local_size, scale)
  for (int i = 0; i < local_size; i++) {
    output[i] = static_cast<unsigned char>(static_cast<float>(local_vec[i] - min) * scale);
  }

  MPI_Allgatherv(output.data(), local_size, MPI_UNSIGNED_CHAR, GetOutput().data(), sendcounts.data(), offsets.data(),
                 MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
  return true;
}

bool LukinITestTaskALL::PostProcessingImpl() {
  return true;
}

}  // namespace lukin_i_ench_contr_lin_hist
