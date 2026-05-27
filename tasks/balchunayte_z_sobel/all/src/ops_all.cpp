#include "balchunayte_z_sobel/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <thread>
#include <vector>

#include "balchunayte_z_sobel/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"

namespace balchunayte_z_sobel {

namespace {

int ConvertPixelToGray(const Pixel &pixel_value) {
  return (77 * static_cast<int>(pixel_value.r) + 150 * static_cast<int>(pixel_value.g) +
          29 * static_cast<int>(pixel_value.b)) >>
         8;
}

void ProcessRows(const Image &input_image, std::vector<int> &output_data, int row_begin, int row_end) {
  const int image_width = input_image.width;
  const auto image_width_size = static_cast<size_t>(image_width);

  for (int row_index = row_begin; row_index < row_end; ++row_index) {
    for (int col_index = 1; col_index < image_width - 1; ++col_index) {
      const size_t index_top_left =
          (static_cast<size_t>(row_index - 1) * image_width_size) + static_cast<size_t>(col_index - 1);
      const size_t index_top_middle =
          (static_cast<size_t>(row_index - 1) * image_width_size) + static_cast<size_t>(col_index);
      const size_t index_top_right =
          (static_cast<size_t>(row_index - 1) * image_width_size) + static_cast<size_t>(col_index + 1);

      const size_t index_middle_left =
          (static_cast<size_t>(row_index) * image_width_size) + static_cast<size_t>(col_index - 1);
      const size_t index_middle_right =
          (static_cast<size_t>(row_index) * image_width_size) + static_cast<size_t>(col_index + 1);

      const size_t index_bottom_left =
          (static_cast<size_t>(row_index + 1) * image_width_size) + static_cast<size_t>(col_index - 1);
      const size_t index_bottom_middle =
          (static_cast<size_t>(row_index + 1) * image_width_size) + static_cast<size_t>(col_index);
      const size_t index_bottom_right =
          (static_cast<size_t>(row_index + 1) * image_width_size) + static_cast<size_t>(col_index + 1);

      const int gray_top_left = ConvertPixelToGray(input_image.data[index_top_left]);
      const int gray_top_middle = ConvertPixelToGray(input_image.data[index_top_middle]);
      const int gray_top_right = ConvertPixelToGray(input_image.data[index_top_right]);

      const int gray_middle_left = ConvertPixelToGray(input_image.data[index_middle_left]);
      const int gray_middle_right = ConvertPixelToGray(input_image.data[index_middle_right]);

      const int gray_bottom_left = ConvertPixelToGray(input_image.data[index_bottom_left]);
      const int gray_bottom_middle = ConvertPixelToGray(input_image.data[index_bottom_middle]);
      const int gray_bottom_right = ConvertPixelToGray(input_image.data[index_bottom_right]);

      const int gradient_x = (-gray_top_left + gray_top_right) + (-2 * gray_middle_left + (2 * gray_middle_right)) +
                             (-gray_bottom_left + gray_bottom_right);

      const int gradient_y = (gray_top_left + (2 * gray_top_middle) + gray_top_right) +
                             (-gray_bottom_left - (2 * gray_bottom_middle) - gray_bottom_right);

      const int magnitude = std::abs(gradient_x) + std::abs(gradient_y);

      const size_t output_index = (static_cast<size_t>(row_index) * image_width_size) + static_cast<size_t>(col_index);
      output_data[output_index] = magnitude;
    }
  }
}

void ProcessRowsOMP(const Image &input_image, std::vector<int> &output_data, int row_begin, int row_end) {
#pragma omp parallel for default(none) shared(input_image, output_data, row_begin, row_end) schedule(static)
  for (int row_index = row_begin; row_index < row_end; ++row_index) {
    ProcessRows(input_image, output_data, row_index, row_index + 1);
  }
}

void ProcessRowsTBB(const Image &input_image, std::vector<int> &output_data, int row_begin, int row_end) {
  oneapi::tbb::parallel_for(row_begin, row_end, [&input_image, &output_data](int row_index) {
    ProcessRows(input_image, output_data, row_index, row_index + 1);
  });
}

void ProcessRowsSTL(const Image &input_image, std::vector<int> &output_data, int row_begin, int row_end) {
  const int local_row_count = row_end - row_begin;

  if (local_row_count <= 0) {
    return;
  }

  const int requested_thread_count = std::max(1, ppc::util::GetNumThreads());
  const int actual_thread_count = std::min(requested_thread_count, local_row_count);

  std::vector<std::thread> worker_threads;
  worker_threads.reserve(static_cast<size_t>(actual_thread_count));

  const int rows_per_thread = local_row_count / actual_thread_count;
  const int remaining_rows = local_row_count % actual_thread_count;

  int start_row_index = row_begin;

  for (int thread_index = 0; thread_index < actual_thread_count; ++thread_index) {
    const int extra_row_count = static_cast<int>(thread_index < remaining_rows);
    const int end_row_index = start_row_index + rows_per_thread + extra_row_count;

    worker_threads.emplace_back(ProcessRows, std::cref(input_image), std::ref(output_data), start_row_index,
                                end_row_index);

    start_row_index = end_row_index;
  }

  for (auto &worker_thread : worker_threads) {
    worker_thread.join();
  }
}

}  // namespace

BalchunayteZSobelOpALL::BalchunayteZSobelOpALL(const InType &input_image) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_image;
  GetOutput().clear();
}

bool BalchunayteZSobelOpALL::ValidationImpl() {
  const auto &input_image = GetInput();

  if (input_image.width <= 0 || input_image.height <= 0) {
    return false;
  }

  const auto expected_size = static_cast<size_t>(input_image.width) * static_cast<size_t>(input_image.height);

  if (input_image.data.size() != expected_size) {
    return false;
  }

  return GetOutput().empty();
}

bool BalchunayteZSobelOpALL::PreProcessingImpl() {
  const auto &input_image = GetInput();
  GetOutput().assign(static_cast<size_t>(input_image.width) * static_cast<size_t>(input_image.height), 0);
  return true;
}

bool BalchunayteZSobelOpALL::RunImpl() {
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const auto &input_image = GetInput();
  auto &output_data = GetOutput();

  if (input_image.width < 3 || input_image.height < 3) {
    MPI_Barrier(MPI_COMM_WORLD);
    return rank >= 0;
  }

  const int first_inner_row = 1;
  const int last_inner_row = input_image.height - 1;
  const int inner_row_count = last_inner_row - first_inner_row;

  const int first_border = first_inner_row + (inner_row_count / 3);
  const int second_border = first_inner_row + ((2 * inner_row_count) / 3);

  ProcessRowsOMP(input_image, output_data, first_inner_row, first_border);
  ProcessRowsSTL(input_image, output_data, first_border, second_border);
  ProcessRowsTBB(input_image, output_data, second_border, last_inner_row);

  MPI_Barrier(MPI_COMM_WORLD);

  return rank >= 0;
}

bool BalchunayteZSobelOpALL::PostProcessingImpl() {
  return true;
}

}  // namespace balchunayte_z_sobel
