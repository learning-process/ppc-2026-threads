#include "fatehov_k_gaussian/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "fatehov_k_gaussian/common/include/common.hpp"

namespace fatehov_k_gaussian {

namespace {

float ConvolvePixel(const std::vector<uint8_t> &data, const std::vector<float> &kernel, int kernel_size, int half,
                    int w, int h, int ch, int y_coord, int x_coord, int c_coord) {
  float res = 0.0F;

  for (int ky = -half; ky <= half; ++ky) {
    for (int kx = -half; kx <= half; ++kx) {
      const int ny = std::clamp(y_coord + ky, 0, h - 1);
      const int nx = std::clamp(x_coord + kx, 0, w - 1);
      const float weight = kernel[((ky + half) * kernel_size) + (kx + half)];
      res += static_cast<float>(data[((ny * w + nx) * ch) + c_coord]) * weight;
    }
  }

  return res;
}

void ProcessPixel(const std::vector<uint8_t> &src_data, std::vector<uint8_t> &dst_data,
                  const std::vector<float> &kernel, int kernel_size, int half, int w, int h, int ch, int y_coord,
                  int x_coord, int row_offset) {
  for (int c_coord = 0; c_coord < ch; ++c_coord) {
    const float res = ConvolvePixel(src_data, kernel, kernel_size, half, w, h, ch, y_coord, x_coord, c_coord);
    dst_data[(((y_coord - row_offset) * w + x_coord) * ch) + c_coord] =
        static_cast<uint8_t>(std::clamp(res, 0.0F, 255.0F));
  }
}

void ProcessRowRangeOMP(const std::vector<uint8_t> &src_data, std::vector<uint8_t> &dst_data,
                        const std::vector<float> &kernel, int kernel_size, int half, int w, int h, int ch,
                        int row_begin, int row_end) {
#pragma omp parallel for default(none) \
    shared(src_data, dst_data, kernel, kernel_size, half, w, h, ch, row_begin, row_end) schedule(static)
  for (int y_coord = row_begin; y_coord < row_end; ++y_coord) {
    for (int x_coord = 0; x_coord < w; ++x_coord) {
      ProcessPixel(src_data, dst_data, kernel, kernel_size, half, w, h, ch, y_coord, x_coord, row_begin);
    }
  }
}

void ComputeDistribution(int h, int size, std::vector<int> &recv_counts, std::vector<int> &displacements, int w,
                         int ch) {
  const int rows_per_proc = h / size;
  const int remainder = h % size;

  for (int i = 0; i < size; ++i) {
    const int r_begin = (i * rows_per_proc) + std::min(i, remainder);
    const int r_end = r_begin + rows_per_proc + (i < remainder ? 1 : 0);
    recv_counts[i] = (r_end - r_begin) * w * ch;
    displacements[i] = r_begin * w * ch;
  }
}

}  // namespace

FatehovKGaussianALL::FatehovKGaussianALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FatehovKGaussianALL::ValidationImpl() {
  const auto &input = GetInput();
  return input.image.width > 0 && input.image.height > 0 && input.image.channels > 0 && !input.image.data.empty() &&
         input.sigma > 0.0F;
}

bool FatehovKGaussianALL::PreProcessingImpl() {
  const auto &input = GetInput();
  const float sigma = input.sigma;

  kernel_size_ = (2 * static_cast<int>(std::ceil(3.0F * sigma))) + 1;
  kernel_.resize(static_cast<std::size_t>(kernel_size_) * kernel_size_);

  const int half = kernel_size_ / 2;
  const float two_sigma_sq = 2.0F * sigma * sigma;
  float sum = 0.0F;

  for (int i = -half; i <= half; ++i) {
    for (int j = -half; j <= half; ++j) {
      const float val = std::exp(-(static_cast<float>((i * i) + (j * j))) / two_sigma_sq);
      kernel_[((i + half) * kernel_size_) + (j + half)] = val;
      sum += val;
    }
  }

  for (float &val : kernel_) {
    val /= sum;
  }

  GetOutput() = Image(input.image.width, input.image.height, input.image.channels);

  return true;
}

bool FatehovKGaussianALL::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const auto &input = GetInput();
  auto &output = GetOutput();
  const int w = static_cast<int>(input.image.width);
  const int h = static_cast<int>(input.image.height);
  const int ch = static_cast<int>(input.image.channels);
  const int half = kernel_size_ / 2;
  const int kernel_size = kernel_size_;
  const auto &kernel = kernel_;

  std::vector<uint8_t> img_data(input.image.data);
  MPI_Bcast(img_data.data(), static_cast<int>(img_data.size()), MPI_UINT8_T, 0, MPI_COMM_WORLD);

  const int rows_per_proc = h / size;
  const int remainder = h % size;
  const int row_begin = (rank * rows_per_proc) + std::min(rank, remainder);
  const int row_end = row_begin + rows_per_proc + (rank < remainder ? 1 : 0);
  const int local_rows = row_end - row_begin;

  std::vector<uint8_t> local_output(static_cast<std::size_t>(local_rows) * w * ch);

  ProcessRowRangeOMP(img_data, local_output, kernel, kernel_size, half, w, h, ch, row_begin, row_end);

  std::vector<int> recv_counts(size);
  std::vector<int> displacements(size);
  ComputeDistribution(h, size, recv_counts, displacements, w, ch);

  MPI_Gatherv(local_output.data(), local_rows * w * ch, MPI_UINT8_T, output.data.data(), recv_counts.data(),
              displacements.data(), MPI_UINT8_T, 0, MPI_COMM_WORLD);

  return true;
}

bool FatehovKGaussianALL::PostProcessingImpl() {
  return true;
}

}  // namespace fatehov_k_gaussian
