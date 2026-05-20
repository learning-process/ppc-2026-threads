#include "kopilov_d_vertical_gauss_filter/all/include/ops_all.hpp"

#include <mpi.h>
#include <oneapi/tbb/blocked_range2d.h>
#include <oneapi/tbb/parallel_for.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"

namespace kopilov_d_vertical_gauss_filter {

namespace {
const int kDivisor = 16;
const std::array<std::array<int, 3>, 3> kGaussKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t GetPixelWithHalo(const uint8_t *local_data, int col_with_halo, int row, int lw_with_halo, int height) {
  int cur_row = row;
  if (cur_row < 0) {
    cur_row = -cur_row - 1;
  } else if (cur_row >= height) {
    cur_row = (2 * height) - cur_row - 1;
  }
  auto idx = (static_cast<size_t>(cur_row) * static_cast<size_t>(lw_with_halo)) + static_cast<size_t>(col_with_halo);
  return local_data[idx];
}

uint8_t ApplyFilter(const uint8_t *data, int col, int row, int lw_halo, int height) {
  int pixel_sum = 0;
  for (size_t ky = 0; ky < 3; ++ky) {
    for (size_t kx = 0; kx < 3; ++kx) {
      int cur_col = col + static_cast<int>(kx);
      int cur_row = row + static_cast<int>(ky) - 1;
      pixel_sum += kGaussKernel.at(ky).at(kx) * GetPixelWithHalo(data, cur_col, cur_row, lw_halo, height);
    }
  }
  return static_cast<uint8_t>(pixel_sum / kDivisor);
}

int FindLeftNeighbor(int rank, const std::vector<int> &counts) {
  for (int p_idx = rank - 1; p_idx >= 0; --p_idx) {
    if (counts.at(static_cast<size_t>(p_idx)) > 0) {
      return p_idx;
    }
  }
  return MPI_PROC_NULL;
}

int FindRightNeighbor(int rank, int size, const std::vector<int> &counts) {
  for (int p_idx = rank + 1; p_idx < size; ++p_idx) {
    if (counts.at(static_cast<size_t>(p_idx)) > 0) {
      return p_idx;
    }
  }
  return MPI_PROC_NULL;
}

void ExchangeHalo(int world_rank, int world_size, int height, int local_w, int lw_with_halo,
                  std::vector<uint8_t> &local_input, const std::vector<int> &counts) {
  if (local_w <= 0) {
    return;
  }

  std::vector<uint8_t> send_left(static_cast<size_t>(height));
  std::vector<uint8_t> send_right(static_cast<size_t>(height));
  std::vector<uint8_t> recv_left(static_cast<size_t>(height));
  std::vector<uint8_t> recv_right(static_cast<size_t>(height));

  for (int row_idx = 0; row_idx < height; ++row_idx) {
    send_left[static_cast<size_t>(row_idx)] =
        local_input[(static_cast<size_t>(row_idx) * static_cast<size_t>(lw_with_halo)) + 1];
    send_right[static_cast<size_t>(row_idx)] =
        local_input[(static_cast<size_t>(row_idx) * static_cast<size_t>(lw_with_halo)) + static_cast<size_t>(local_w)];
  }

  int left_proc = FindLeftNeighbor(world_rank, counts);
  int right_proc = FindRightNeighbor(world_rank, world_size, counts);

  MPI_Sendrecv(send_left.data(), height, MPI_UNSIGNED_CHAR, left_proc, 1, recv_right.data(), height, MPI_UNSIGNED_CHAR,
               right_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Sendrecv(send_right.data(), height, MPI_UNSIGNED_CHAR, right_proc, 2, recv_left.data(), height, MPI_UNSIGNED_CHAR,
               left_proc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  for (int row_idx = 0; row_idx < height; ++row_idx) {
    auto base_idx = static_cast<size_t>(row_idx) * static_cast<size_t>(lw_with_halo);
    if (left_proc != MPI_PROC_NULL) {
      local_input[base_idx + 0] = recv_left[static_cast<size_t>(row_idx)];
    } else {
      local_input[base_idx + 0] = local_input[base_idx + 1];
    }

    if (right_proc != MPI_PROC_NULL) {
      local_input[base_idx + static_cast<size_t>(local_w + 1)] = recv_right[static_cast<size_t>(row_idx)];
    } else {
      local_input[base_idx + static_cast<size_t>(local_w + 1)] = local_input[base_idx + static_cast<size_t>(local_w)];
    }
  }
}

void MasterDistribute(int world_size, int width, int height, const std::vector<int> &counts,
                      const std::vector<int> &displs, const uint8_t *full_data, uint8_t *local_data, int lw_with_halo) {
  int p0_w = counts.at(0);
  for (int row_idx = 0; row_idx < height; ++row_idx) {
    for (int col_idx = 0; col_idx < p0_w; ++col_idx) {
      local_data[(static_cast<size_t>(row_idx) * static_cast<size_t>(lw_with_halo)) +
                 static_cast<size_t>(col_idx + 1)] =
          full_data[(static_cast<size_t>(row_idx) * static_cast<size_t>(width)) + static_cast<size_t>(col_idx)];
    }
  }

  for (int proc_idx = 1; proc_idx < world_size; ++proc_idx) {
    int cur_w = counts.at(static_cast<size_t>(proc_idx));
    if (cur_w == 0) {
      continue;
    }
    std::vector<uint8_t> buf(static_cast<size_t>(cur_w) * static_cast<size_t>(height));
    for (int row_idx = 0; row_idx < height; ++row_idx) {
      for (int col_idx = 0; col_idx < cur_w; ++col_idx) {
        buf[(static_cast<size_t>(row_idx) * static_cast<size_t>(cur_w)) + static_cast<size_t>(col_idx)] =
            full_data[(static_cast<size_t>(row_idx) * static_cast<size_t>(width)) +
                      static_cast<size_t>(displs.at(static_cast<size_t>(proc_idx)) + col_idx)];
      }
    }
    MPI_Send(buf.data(), static_cast<int>(buf.size()), MPI_UNSIGNED_CHAR, proc_idx, 0, MPI_COMM_WORLD);
  }
}

void WorkerReceive(int my_w, int height, uint8_t *local_data, int lw_with_halo) {
  if (my_w <= 0) {
    return;
  }
  std::vector<uint8_t> buf(static_cast<size_t>(my_w) * static_cast<size_t>(height));
  MPI_Recv(buf.data(), static_cast<int>(buf.size()), MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  for (int row_idx = 0; row_idx < height; ++row_idx) {
    for (int col_idx = 0; col_idx < my_w; ++col_idx) {
      local_data[(static_cast<size_t>(row_idx) * static_cast<size_t>(lw_with_halo)) +
                 static_cast<size_t>(col_idx + 1)] =
          buf[(static_cast<size_t>(row_idx) * static_cast<size_t>(my_w)) + static_cast<size_t>(col_idx)];
    }
  }
}

void MasterGather(int world_size, int width, int height, const std::vector<int> &counts, const std::vector<int> &displs,
                  const std::vector<uint8_t> &local_out, uint8_t *full_data) {
  int p0_w = counts.at(0);
  for (int row_idx = 0; row_idx < height; ++row_idx) {
    for (int col_idx = 0; col_idx < p0_w; ++col_idx) {
      full_data[(static_cast<size_t>(row_idx) * static_cast<size_t>(width)) + static_cast<size_t>(col_idx)] =
          local_out[(static_cast<size_t>(row_idx) * static_cast<size_t>(p0_w)) + static_cast<size_t>(col_idx)];
    }
  }

  for (int proc_idx = 1; proc_idx < world_size; ++proc_idx) {
    int cur_w = counts.at(static_cast<size_t>(proc_idx));
    if (cur_w == 0) {
      continue;
    }
    std::vector<uint8_t> buf(static_cast<size_t>(cur_w) * static_cast<size_t>(height));
    MPI_Recv(buf.data(), static_cast<int>(buf.size()), MPI_UNSIGNED_CHAR, proc_idx, 3, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (int row_idx = 0; row_idx < height; ++row_idx) {
      for (int col_idx = 0; col_idx < cur_w; ++col_idx) {
        full_data[(static_cast<size_t>(row_idx) * static_cast<size_t>(width)) +
                  static_cast<size_t>(displs.at(static_cast<size_t>(proc_idx)) + col_idx)] =
            buf[(static_cast<size_t>(row_idx) * static_cast<size_t>(cur_w)) + static_cast<size_t>(col_idx)];
      }
    }
  }
}

void WorkerGather(int local_w, const std::vector<uint8_t> &local_out) {
  if (local_w <= 0) {
    return;
  }
  MPI_Send(local_out.data(), static_cast<int>(local_out.size()), MPI_UNSIGNED_CHAR, 0, 3, MPI_COMM_WORLD);
}

void ComputeTBB(int local_w, int global_height, int lw_halo, const uint8_t *in_data, uint8_t *out_data) {
  if (local_w <= 0) {
    return;
  }
  tbb::parallel_for(tbb::blocked_range2d<int>(0, global_height, 0, local_w), [&](const tbb::blocked_range2d<int> &r) {
    for (int row_idx = r.rows().begin(); row_idx != r.rows().end(); ++row_idx) {
      for (int col_idx = r.cols().begin(); col_idx != r.cols().end(); ++col_idx) {
        out_data[(static_cast<size_t>(row_idx) * static_cast<size_t>(local_w)) + static_cast<size_t>(col_idx)] =
            ApplyFilter(in_data, col_idx, row_idx, lw_halo, global_height);
      }
    }
  });
}

}  // namespace

KopilovDVerticalGaussFilterALL::KopilovDVerticalGaussFilterALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool KopilovDVerticalGaussFilterALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int res = 0;
  if (rank == 0) {
    const auto &in = GetInput();
    if (in.width > 0 && in.height > 0 &&
        in.data.size() == static_cast<size_t>(in.width) * static_cast<size_t>(in.height)) {
      res = 1;
    }
  }
  MPI_Bcast(&res, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return res == 1;
}

bool KopilovDVerticalGaussFilterALL::PreProcessingImpl() {
  return true;
}

bool KopilovDVerticalGaussFilterALL::RunImpl() {
  int world_size = 0;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int global_width = 0;
  int global_height = 0;
  if (world_rank == 0) {
    global_width = GetInput().width;
    global_height = GetInput().height;
  }
  MPI_Bcast(&global_width, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_height, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (global_width <= 0 || global_height <= 0) {
    return true;
  }

  std::vector<int> counts(static_cast<size_t>(world_size));
  std::vector<int> displs(static_cast<size_t>(world_size));
  int offset = 0;
  for (int i = 0; i < world_size; ++i) {
    counts[static_cast<size_t>(i)] = (global_width / world_size) + (i < (global_width % world_size) ? 1 : 0);
    displs[static_cast<size_t>(i)] = offset;
    offset += counts[static_cast<size_t>(i)];
  }

  int local_w = counts.at(static_cast<size_t>(world_rank));
  int lw_halo = local_w + 2;
  std::vector<uint8_t> l_in(static_cast<size_t>(lw_halo) * static_cast<size_t>(global_height), 0);
  std::vector<uint8_t> l_out(static_cast<size_t>(local_w) * static_cast<size_t>(global_height));

  if (world_rank == 0) {
    MasterDistribute(world_size, global_width, global_height, counts, displs, GetInput().data.data(), l_in.data(),
                     lw_halo);
  } else {
    WorkerReceive(local_w, global_height, l_in.data(), lw_halo);
  }

  ExchangeHalo(world_rank, world_size, global_height, local_w, lw_halo, l_in, counts);
  ComputeTBB(local_w, global_height, lw_halo, l_in.data(), l_out.data());

  if (world_rank == 0) {
    GetOutput().width = global_width;
    GetOutput().height = global_height;
    GetOutput().data.assign(static_cast<size_t>(global_width) * static_cast<size_t>(global_height), 0);
    MasterGather(world_size, global_width, global_height, counts, displs, l_out, GetOutput().data.data());
  } else {
    WorkerGather(local_w, l_out);
  }

  return true;
}

bool KopilovDVerticalGaussFilterALL::PostProcessingImpl() {
  return true;
}

}  // namespace kopilov_d_vertical_gauss_filter
