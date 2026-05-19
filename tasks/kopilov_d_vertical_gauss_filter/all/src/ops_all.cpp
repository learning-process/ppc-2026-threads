#include "kopilov_d_vertical_gauss_filter/all/include/ops_all.hpp"

#include <mpi.h>
#include <oneapi/tbb/blocked_range2d.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
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
  int cur_row = std::clamp(row, 0, height - 1);
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

void ExchangeHalo(int world_rank, int world_size, int height, int local_w, int lw_with_halo,
                  std::vector<uint8_t> &local_input) {
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

  int left_proc = (world_rank > 0) ? world_rank - 1 : MPI_PROC_NULL;
  int right_proc = (world_rank < world_size - 1) ? world_rank + 1 : MPI_PROC_NULL;

  MPI_Sendrecv(send_left.data(), height, MPI_UNSIGNED_CHAR, left_proc, 1, recv_right.data(), height, MPI_UNSIGNED_CHAR,
               right_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Sendrecv(send_right.data(), height, MPI_UNSIGNED_CHAR, right_proc, 2, recv_left.data(), height, MPI_UNSIGNED_CHAR,
               left_proc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  for (int row_idx = 0; row_idx < height; ++row_idx) {
    auto base_idx = static_cast<size_t>(row_idx) * static_cast<size_t>(lw_with_halo);
    local_input[base_idx + 0] =
        (left_proc != MPI_PROC_NULL) ? recv_left[static_cast<size_t>(row_idx)] : local_input[base_idx + 1];
    local_input[base_idx + static_cast<size_t>(local_w + 1)] =
        (right_proc != MPI_PROC_NULL) ? recv_right[static_cast<size_t>(row_idx)]
                                      : local_input[base_idx + static_cast<size_t>(local_w)];
  }
}

void MasterDistribute(int world_size, int width, int height, const std::vector<int> &counts,
                      const std::vector<int> &displs, const uint8_t *full_data, uint8_t *local_data, int lw_with_halo) {
  int p0_w = counts.at(0);
  for (int r_idx = 0; r_idx < height; ++r_idx) {
    for (int c_idx = 0; c_idx < p0_w; ++c_idx) {
      local_data[(static_cast<size_t>(r_idx) * static_cast<size_t>(lw_with_halo)) + static_cast<size_t>(c_idx + 1)] =
          full_data[(static_cast<size_t>(r_idx) * static_cast<size_t>(width)) + static_cast<size_t>(c_idx)];
    }
  }

  for (int p_idx = 1; p_idx < world_size; ++p_idx) {
    int p_w = counts.at(static_cast<size_t>(p_idx));
    std::vector<uint8_t> buf(static_cast<size_t>(p_w) * static_cast<size_t>(height));
    for (int r_idx = 0; r_idx < height; ++r_idx) {
      for (int c_idx = 0; c_idx < p_w; ++c_idx) {
        buf[(static_cast<size_t>(r_idx) * static_cast<size_t>(p_w)) + static_cast<size_t>(c_idx)] =
            full_data[(static_cast<size_t>(r_idx) * static_cast<size_t>(width)) +
                      static_cast<size_t>(displs.at(static_cast<size_t>(p_idx)) + c_idx)];
      }
    }
    MPI_Send(buf.data(), static_cast<int>(buf.size()), MPI_UNSIGNED_CHAR, p_idx, 0, MPI_COMM_WORLD);
  }
}

void WorkerReceive(int my_w, int height, uint8_t *local_data, int lw_with_halo) {
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
  for (int r_idx = 0; r_idx < height; ++r_idx) {
    for (int c_idx = 0; c_idx < p0_w; ++c_idx) {
      full_data[(static_cast<size_t>(r_idx) * static_cast<size_t>(width)) + static_cast<size_t>(c_idx)] =
          local_out[(static_cast<size_t>(r_idx) * static_cast<size_t>(p0_w)) + static_cast<size_t>(c_idx)];
    }
  }

  for (int p_idx = 1; p_idx < world_size; ++p_idx) {
    int p_w = counts.at(static_cast<size_t>(p_idx));
    std::vector<uint8_t> buf(static_cast<size_t>(p_w) * static_cast<size_t>(height));
    MPI_Recv(buf.data(), static_cast<int>(buf.size()), MPI_UNSIGNED_CHAR, p_idx, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int r_idx = 0; r_idx < height; ++r_idx) {
      for (int c_idx = 0; c_idx < p_w; ++c_idx) {
        full_data[(static_cast<size_t>(r_idx) * static_cast<size_t>(width)) +
                  static_cast<size_t>(displs.at(static_cast<size_t>(p_idx)) + c_idx)] =
            buf[(static_cast<size_t>(r_idx) * static_cast<size_t>(p_w)) + static_cast<size_t>(c_idx)];
      }
    }
  }
}
}  // namespace

KopilovDVerticalGaussFilterALL::KopilovDVerticalGaussFilterALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool KopilovDVerticalGaussFilterALL::ValidationImpl() {
  int rank_v = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_v);
  if (rank_v == 0) {
    const auto &input_ref = GetInput();
    return input_ref.width > 0 && input_ref.height > 0 &&
           input_ref.data.size() == static_cast<size_t>(input_ref.width) * static_cast<size_t>(input_ref.height);
  }
  return true;
}

bool KopilovDVerticalGaussFilterALL::PreProcessingImpl() {
  return true;
}

bool KopilovDVerticalGaussFilterALL::RunImpl() {
  int world_size = 0;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int g_w = 0;
  int g_h = 0;
  if (world_rank == 0) {
    g_w = GetInput().width;
    g_h = GetInput().height;
  }
  MPI_Bcast(&g_w, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (g_w <= 0 || g_h <= 0) {
    return true;
  }
  std::vector<int> counts(static_cast<size_t>(world_size));
  std::vector<int> displs(static_cast<size_t>(world_size));
  int offset_val = 0;
  for (int idx = 0; idx < world_size; ++idx) {
    counts[static_cast<size_t>(idx)] = (g_w / world_size) + (idx < (g_w % world_size) ? 1 : 0);
    displs[static_cast<size_t>(idx)] = offset_val;
    offset_val += counts[static_cast<size_t>(idx)];
  }
  int local_w = counts.at(static_cast<size_t>(world_rank));
  int lw_halo = local_w + 2;
  std::vector<uint8_t> l_in(static_cast<size_t>(lw_halo) * static_cast<size_t>(g_h), 0);
  std::vector<uint8_t> l_out(static_cast<size_t>(local_w) * static_cast<size_t>(g_h));

  if (world_rank == 0) {
    MasterDistribute(world_size, g_w, g_h, counts, displs, GetInput().data.data(), l_in.data(), lw_halo);
  } else {
    WorkerReceive(local_w, g_h, l_in.data(), lw_halo);
  }
  ExchangeHalo(world_rank, world_size, g_h, local_w, lw_halo, l_in);
  tbb::parallel_for(tbb::blocked_range2d<int>(0, g_h, 0, local_w), [&](const tbb::blocked_range2d<int> &r) {
    for (int r_idx = r.rows().begin(); r_idx != r.rows().end(); ++r_idx) {
      for (int c_idx = r.cols().begin(); c_idx != r.cols().end(); ++c_idx) {
        l_out[(static_cast<size_t>(r_idx) * static_cast<size_t>(local_w)) + static_cast<size_t>(c_idx)] =
            ApplyFilter(l_in.data(), c_idx, r_idx, lw_halo, g_h);
      }
    }
  });
  if (world_rank == 0) {
    GetOutput().width = g_w;
    GetOutput().height = g_h;
    GetOutput().data.resize(static_cast<size_t>(g_w) * static_cast<size_t>(g_h));
    MasterGather(world_size, g_w, g_h, counts, displs, l_out, GetOutput().data.data());
  } else {
    MPI_Send(l_out.data(), static_cast<int>(l_out.size()), MPI_UNSIGNED_CHAR, 0, 3, MPI_COMM_WORLD);
  }
  return true;
}

bool KopilovDVerticalGaussFilterALL::PostProcessingImpl() {
  return true;
}

}  // namespace kopilov_d_vertical_gauss_filter
