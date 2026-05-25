#include "fedoseev_linear_image_filtering_vertical/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"
#include "util/include/util.hpp"

namespace fedoseev_linear_image_filtering_vertical {

namespace {
int ComputePixel(const std::vector<int> &img, int w, int h, int row, int col,
                 const std::array<std::array<int, 3>, 3> &kernel) {
  const int kernel_sum = 16;
  int sum = 0;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      int px = col + kx;
      int py = row + ky;
      px = std::clamp(px, 0, w - 1);
      py = std::clamp(py, 0, h - 1);
      sum += img[(static_cast<size_t>(py) * static_cast<size_t>(w)) + static_cast<size_t>(px)] *
             kernel.at(ky + 1).at(kx + 1);
    }
  }
  return sum / kernel_sum;
}
}  // namespace

LinearImageFilteringVerticalAll::LinearImageFilteringVerticalAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = InType{};
}

bool LinearImageFilteringVerticalAll::ValidationImpl() {
  const InType &input = GetInput();
  if (input.width < 3 || input.height < 3) {
    return false;
  }
  if (input.data.size() != static_cast<size_t>(input.width) * static_cast<size_t>(input.height)) {
    return false;
  }
  int initialized = 0;
  MPI_Initialized(&initialized);
  return initialized != 0;
}

bool LinearImageFilteringVerticalAll::PreProcessingImpl() {
  return true;
}

bool LinearImageFilteringVerticalAll::RunImpl() {
  int rank = 0;
  int size = 1;
  bool is_mpi = ppc::util::IsUnderMpirun();
  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  int w = 0;
  int h = 0;
  std::vector<int> img;
  if (rank == 0) {
    const InType &input = GetInput();
    w = input.width;
    h = input.height;
    img = input.data;
  }

  if (is_mpi) {
    MPI_Bcast(&w, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
      img.resize(static_cast<size_t>(w) * static_cast<size_t>(h));
    }
    MPI_Bcast(img.data(), static_cast<int>(img.size()), MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    const InType &input = GetInput();
    w = input.width;
    h = input.height;
    img = input.data;
  }

  const std::array<std::array<int, 3>, 3> kernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
  std::vector<int> result(static_cast<size_t>(w) * static_cast<size_t>(h), 0);

  for (int row = 0; row < h; ++row) {
    for (int col = 0; col < w; ++col) {
      result[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)] =
          ComputePixel(img, w, h, row, col, kernel);
    }
  }

  OutType out;
  out.width = w;
  out.height = h;
  out.data = std::move(result);
  GetOutput() = out;

  return true;
}

bool LinearImageFilteringVerticalAll::PostProcessingImpl() {
  return true;
}

}  // namespace fedoseev_linear_image_filtering_vertical
