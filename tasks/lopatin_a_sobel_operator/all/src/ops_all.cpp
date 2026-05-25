#include "lopatin_a_sobel_operator/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "lopatin_a_sobel_operator/common/include/common.hpp"
#include "util/include/util.hpp"

namespace lopatin_a_sobel_operator {

const std::array<std::array<int, 3>, 3> kSobelX = {std::array<int, 3>{-1, 0, 1}, std::array<int, 3>{-2, 0, 2},
                                                   std::array<int, 3>{-1, 0, 1}};

const std::array<std::array<int, 3>, 3> kSobelY = {std::array<int, 3>{-1, -2, -1}, std::array<int, 3>{0, 0, 0},
                                                   std::array<int, 3>{1, 2, 1}};

void LopatinASobelOperatorALL::RunSobel(const Image &img, std::size_t start, std::size_t end,
                                        lopatin_a_sobel_operator::OutType &output, std::size_t shift) {
  const auto &input_data = img.pixels;

  int proc_num = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

#pragma omp parallel for num_threads(ppc::util::GetNumThreads() / proc_num) schedule(static) default(none) \
    shared(kSobelX, kSobelY, img, input_data, output, start, end, shift)
  for (std::size_t j = start; j < end; ++j) {  // processing only pixels with a full 3 x 3 neighborhood size
    std::size_t out_row = j - shift;
    for (std::size_t i = 1; i < img.width - 1; ++i) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          std::uint8_t pixel = input_data[((j + ky) * img.width) + (i + kx)];
          gx += pixel * kSobelX.at(ky + 1).at(kx + 1);
          gy += pixel * kSobelY.at(ky + 1).at(kx + 1);
        }
      }

      auto magnitude = static_cast<int>(round(std::sqrt((gx * gx) + (gy * gy))));
      output[(out_row * img.width) + i] = (magnitude > img.threshold) ? magnitude : 0;
    }
  }
}

LopatinASobelOperatorALL::LopatinASobelOperatorALL(const InType &in) : h_(in.height), w_(in.width) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool LopatinASobelOperatorALL::ValidationImpl() {
  const auto &input = GetInput();
  return h_ * w_ == input.pixels.size();
}

bool LopatinASobelOperatorALL::PreProcessingImpl() {
  GetOutput().resize(h_ * w_);
  return true;
}

bool LopatinASobelOperatorALL::RunImpl() {
  int proc_num = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
  int proc_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

  const auto &input = GetInput();

  const auto rows = input.height - 2;
  const auto n_procs = std::min(static_cast<std::size_t>(proc_num), rows);

  const std::size_t portion = rows / n_procs;

  lopatin_a_sobel_operator::OutType local_output(portion * input.width, 0);

  if (std::cmp_less(proc_rank, n_procs)) {
    std::size_t start = 1 + (proc_rank * portion);
    std::size_t end = 1 + ((proc_rank + 1) * portion);
    RunSobel(input, start, end, local_output, start);
  }

  auto &output = GetOutput();

  std::vector<int> recvcounts(proc_num, 0);
  std::vector<int> displs(proc_num, 0);

  for (size_t i = 0; i < n_procs; ++i) {
    recvcounts[i] = static_cast<int>(portion * input.width);
    displs[i] = static_cast<int>((1 + i * portion) * input.width);
  }

  int sendcount = std::cmp_less(proc_rank, n_procs) ? static_cast<int>(local_output.size()) : 0;

  MPI_Gatherv(local_output.data(), sendcount, MPI_INT, output.data(), recvcounts.data(), displs.data(), MPI_INT, 0,
              MPI_COMM_WORLD);

  if (proc_rank == 0) {
    const std::size_t tail = rows % n_procs;

    if (tail != 0) {
      std::size_t start = input.height - 1 - tail;
      RunSobel(input, start, input.height - 1, output, 0);
    }
  }

  MPI_Bcast(output.data(), static_cast<int>(output.size()), MPI_INT, 0, MPI_COMM_WORLD);

  return true;
}

bool LopatinASobelOperatorALL::PostProcessingImpl() {
  return true;
}

}  // namespace lopatin_a_sobel_operator
