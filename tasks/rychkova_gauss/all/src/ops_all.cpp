#include "rychkova_gauss/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "rychkova_gauss/common/include/common.hpp"
#include "util/include/util.hpp"

namespace rychkova_gauss {

namespace {
int Mirror(int x, int xmin, int xmax) {
  if (x < xmin) {
    return 1;
  }
  if (x >= xmax) {
    return xmax - 1;
  }
  return x;
};

Pixel ComputePixel(const Image &image, std::size_t x, std::size_t y, std::size_t width, std::size_t height) {
  Pixel result = {.R = 0, .G = 0, .B = 0};
  for (int shift_x = -1; shift_x < 2; shift_x++) {
    for (int shift_y = -1; shift_y < 2; shift_y++) {
      int xn = Mirror(static_cast<int>(x) + shift_x, 0, static_cast<int>(width));
      int yn = Mirror(static_cast<int>(y) + shift_y, 0, static_cast<int>(height));
      auto current = image[yn][xn];
      result.R += static_cast<uint8_t>(static_cast<double>(current.R) * kKernel[shift_x + 1][shift_y + 1]);
      result.G += static_cast<uint8_t>(static_cast<double>(current.G) * kKernel[shift_x + 1][shift_y + 1]);
      result.B += static_cast<uint8_t>(static_cast<double>(current.B) * kKernel[shift_x + 1][shift_y + 1]);
    }
  }
  return result;
}
}  // namespace

RychkovaGaussALL::RychkovaGaussALL(const InType &in) : loutput_({}) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};  // cписок инициализации - пустой вектор - каждый вложенный вектор и пиксели внутри по умолчанию
}

bool RychkovaGaussALL::ValidationImpl() {
  if (GetInput().empty()) {
    return false;
  }
  const auto len = GetInput()[0].size();
  return std::ranges::all_of(GetInput(), [len](const auto &row) { return row.size() == len; });
}

bool RychkovaGaussALL::PreProcessingImpl() {
  const auto &image = GetInput();  // сохран.изоб.
  const auto width = image[0].size();
  const auto height = image.size();
  loutput_.resize(width * height * 3, 0);
  return true;
}

bool RychkovaGaussALL::RunImpl() {
  const auto &image = GetInput();  // сохран.изоб.
  const auto width = image[0].size();
  const auto height = image.size();  // сайз от имаге хранит количествo строк
  GetOutput() = Image(height, std::vector<Pixel>(width, Pixel(0, 0, 0)));
  int n = 0;
  int idx = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  MPI_Comm_rank(MPI_COMM_WORLD, &idx);
  size_t row_per_thread = height / n;
  size_t remainder = height % n;
  std::vector<int> counts(n, 0);
  std::vector<int> displs(n, 0);
  for (int i = 0; i < n; i++) {
    size_t start = (i * row_per_thread) + std::min(static_cast<int>(remainder), i);
    size_t end = ((i + 1) * row_per_thread) + std::min(static_cast<int>(remainder), i + 1);
    counts[i] = static_cast<int>((end - start) * 3 * width);
    displs[i] = static_cast<int>(start * 3 * width);
  }
  size_t start = displs[idx];
  size_t end = start + counts[idx];
#pragma omp parallel for shared(width, height, image, start, end) default(none) num_threads(ppc::util::GetNumThreads())
  for (std::size_t j = start; j < end; j++) {
    size_t x = (j / 3) % width;
    size_t y = (j / 3) / width;
    auto px = ComputePixel(image, x, y, width, height);
    loutput_[j] = px.R;
    loutput_[j + 1] = px.G;
    loutput_[j + 2] = px.B;
  }
  MPI_Gatherv(loutput_.data(), counts[idx], MPI_UINT8_T, loutput_.data(), counts.data(), displs.data(), MPI_UINT8_T, 0,
              MPI_COMM_WORLD);
  return true;
}

bool RychkovaGaussALL::PostProcessingImpl() {
  const auto &image = GetInput();  // сохран.изоб.
  const auto width = image[0].size();
  const auto height = image.size();
  auto &output = GetOutput();
  for (std::size_t j = 0; j < height; j++) {
    for (std::size_t i = 0; i < width; i++) {
      Pixel px(loutput_[((j * width) + 1) * 3], loutput_[(((j * width) + 1) * 3) + 1],
               loutput_[(((j * width) + 1) * 3) + 2]);
      output[j][i] = px;
    }
  }
  return true;
}

}  // namespace rychkova_gauss
