#include "litvyakov_d_shell_sort/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "litvyakov_d_shell_sort/common/include/common.hpp"

namespace litvyakov_d_shell_sort {

void LitvyakovDShellSortALL::BaseShellSort(std::vector<int>::iterator first, std::vector<int>::iterator last) {
  for (std::ptrdiff_t dist = (last - first) / 2; dist > 0; dist /= 2) {
    for (auto i = first + dist; i != last; ++i) {
      for (auto j = i; j - first >= dist && (*j < *(j - dist)); j -= dist) {
        std::swap(*j, *(j - dist));
      }
    }
  }
}

std::vector<std::size_t> LitvyakovDShellSortALL::GetBounds(std::size_t n, std::size_t parts) {
  parts = std::max<std::size_t>(1, std::min(parts, n));

  std::vector<std::size_t> bounds;
  bounds.reserve(parts + 1);
  bounds.push_back(0);

  const std::size_t base = n / parts;
  const std::size_t rem = n % parts;

  for (std::size_t i = 0; i < parts; ++i) {
    bounds.push_back(bounds.back() + base);
    if (i < rem) {
      bounds[i + 1]++;
    }
  }

  return bounds;
}

LitvyakovDShellSortALL::LitvyakovDShellSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    GetInput() = in;
  } else {
    GetInput() = InType();
  }
  GetOutput() = std::vector<int>();
}

bool LitvyakovDShellSortALL::ValidationImpl() {
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    const InType &vec = GetInput();
    return !vec.empty();
  }
  return true;
}

bool LitvyakovDShellSortALL::PreProcessingImpl() {
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    GetOutput() = GetInput();
  }
  return true;
}

bool LitvyakovDShellSortALL::RunImpl() {
  int world_rank = 0;
  int world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  std::vector<int> &vec = GetOutput();

  int n = 0;
  if (world_rank == 0) {
    n = static_cast<int>(vec.size());
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (n <= 1) {
    return true;
  }

  std::vector<int> sendcounts(world_size, 0);
  std::vector<int> displs(world_size, 0);
  int base = n / world_size;
  int rem = n % world_size;
  int current_displ = 0;

  for (int i = 0; i < world_size; ++i) {
    sendcounts[i] = base + (i < rem ? 1 : 0);
    displs[i] = current_displ;
    current_displ += sendcounts[i];
  }

  std::vector<int> local_vec(sendcounts[world_rank]);

  MPI_Scatterv(world_rank == 0 ? vec.data() : nullptr, sendcounts.data(), displs.data(), MPI_INT, local_vec.data(),
               sendcounts[world_rank], MPI_INT, 0, MPI_COMM_WORLD);

  if (!local_vec.empty()) {
    const std::size_t threads = std::max(1, omp_get_max_threads());
    const std::size_t parts_count = std::min<std::size_t>(threads, local_vec.size());
    const auto bounds = GetBounds(local_vec.size(), parts_count);
    int parts_count_t = static_cast<int>(parts_count);

#pragma omp parallel for default(none) shared(local_vec, bounds, parts_count_t) schedule(static)
    for (int i = 0; i < parts_count_t; ++i) {
      const std::size_t l = bounds[i];
      const std::size_t r = bounds[i + 1];
      BaseShellSort(local_vec.begin() + static_cast<std::ptrdiff_t>(l),
                    local_vec.begin() + static_cast<std::ptrdiff_t>(r));
    }

    // cppreference.com:
    //  void inplace_merge( BidirIt first, BidirIt middle, BidirIt last ),
    //  Merges two consecutive sorted ranges [first, middle) and [middle, last) into one sorted range [first, last).
    for (std::size_t i = 1; i < parts_count; ++i) {
      std::inplace_merge(local_vec.begin(), local_vec.begin() + static_cast<std::ptrdiff_t>(bounds[i]),
                         local_vec.begin() + static_cast<std::ptrdiff_t>(bounds[i + 1]));
    }
  }

  MPI_Gatherv(local_vec.data(), sendcounts[world_rank], MPI_INT, world_rank == 0 ? vec.data() : nullptr,
              sendcounts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (world_rank == 0) {
    for (int i = 1; i < world_size; ++i) {
      if (sendcounts[i] > 0) {
        std::inplace_merge(vec.begin(), vec.begin() + static_cast<std::ptrdiff_t>(displs[i]),
                           vec.begin() + static_cast<std::ptrdiff_t>(displs[i] + sendcounts[i]));
      }
    }
  }

  return true;
}

bool LitvyakovDShellSortALL::PostProcessingImpl() {
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  std::vector<int> &vec = GetOutput();
  int n = 0;

  if (world_rank == 0) {
    n = static_cast<int>(vec.size());
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (world_rank != 0) {
    vec.resize(n);
  }

  if (n > 0) {
    MPI_Bcast(vec.data(), n, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return true;
}

}  // namespace litvyakov_d_shell_sort
