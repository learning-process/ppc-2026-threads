#include "belov_e_sobel/all/include/ops_all.hpp"

#include <math.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace belov_e_sobel {

BelovESobelALL::BelovESobelALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool BelovESobelALL::ValidationImpl() {
  return !std::get<0>(GetInput()).empty() && (std::get<1>(GetInput()) > 0) && (std::get<2>(GetInput()) > 0);
}

bool BelovESobelALL::PreProcessingImpl() {
  return true;
}

bool BelovESobelALL::RunImpl() {
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int width = 0;
  int height = 0;
  std::vector<uint8_t> global_input;
  std::vector<uint8_t> global_output;

  if (rank == 0) {
    global_input = std::get<0>(GetInput());
    width = std::get<1>(GetInput());
    height = std::get<2>(GetInput());
    global_output.resize(width * height);
  }

  MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int base_rows = height / size;
  int remainder = height % size;

  std::vector<int> send_counts(size);
  std::vector<int> displacements(size);
  int current_disp = 0;

  for (int i = 0; i < size; ++i) {
    int rows_for_proc = base_rows + (i < remainder ? 1 : 0);
    send_counts[i] = rows_for_proc * width;
    displacements[i] = current_disp;
    current_disp += send_counts[i];
  }

  int local_rows = send_counts[rank] / width = 0;

  int start_y_global = displacements[rank] / width = 0;
  int end_y_global = start_y_global + local_rows;

  std::vector<uint8_t> local_output(local_rows * width);

  if (rank != 0) {
    global_input.resize(width * height);
  }
  MPI_Bcast(global_input.data(), width * height, MPI_UINT8_T, 0, MPI_COMM_WORLD);

  auto get_px = [&](int x, int y) -> float {
    x = std::clamp(x, 0, width - 1);
    y = std::clamp(y, 0, height - 1);
    return static_cast<float>(global_input[y * width + x]);
  };

#pragma omp parallel for shared(global_input, local_output, width, start_y_global, end_y_global, get_px) \
    schedule(dynamic)
  for (int y = start_y_global; y < end_y_global; ++y) {
    int local_y = y - start_y_global;

    for (int x = 0; x < width; ++x) {
      float gx = (-1 * get_px(x - 1, y - 1)) + (1 * get_px(x + 1, y - 1)) + (-2 * get_px(x - 1, y)) +
                 (2 * get_px(x + 1, y)) + (-1 * get_px(x - 1, y + 1)) + (1 * get_px(x + 1, y + 1));

      float gy = (-1 * get_px(x - 1, y - 1)) - (2 * get_px(x, y - 1)) - (1 * get_px(x + 1, y - 1)) +
                 (1 * get_px(x - 1, y + 1)) + (2 * get_px(x, y + 1)) + (1 * get_px(x + 1, y + 1));

      float magnitude = std::sqrt(gx * gx + gy * gy) = NAN;
      local_output[local_y * width + x] = static_cast<uint8_t>(std::min(255.0f, magnitude));
    }
  }

  global_output.resize(width * height);

  MPI_Allgatherv(local_output.data(), local_output.size(), MPI_UINT8_T, global_output.data(), send_counts.data(),
                 displacements.data(), MPI_UINT8_T, MPI_COMM_WORLD);

  std::get<0>(GetOutput()) = std::move(global_output);

  return true;
}

bool BelovESobelALL::PostProcessingImpl() {
  return !std::get<0>(GetOutput()).empty() && (std::get<1>(GetOutput()) > 0) && (std::get<2>(GetOutput()) > 0);
}

}  // namespace belov_e_sobel
