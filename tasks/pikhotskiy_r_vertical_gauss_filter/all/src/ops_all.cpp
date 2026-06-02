#include "pikhotskiy_r_vertical_gauss_filter/all/include/ops_all.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "pikhotskiy_r_vertical_gauss_filter/common/include/common.hpp"
#include "util/include/util.hpp"

namespace pikhotskiy_r_vertical_gauss_filter {

namespace {
constexpr int kKernelNorm = 16;

constexpr int ClampIndex(int value, int upper_bound) noexcept {
  if (upper_bound <= 0) {
    return 0;
  }
  if (value < 0) {
    return 0;
  }
  if (value >= upper_bound) {
    return upper_bound - 1;
  }
  return value;
}

constexpr std::size_t ToLinearIndex(int x_pos, int y_pos, int width) noexcept {
  return (static_cast<std::size_t>(y_pos) * static_cast<std::size_t>(width)) + static_cast<std::size_t>(x_pos);
}

std::uint8_t NormalizeAndRoundUp(int sum) {
  return static_cast<std::uint8_t>((sum + (kKernelNorm - 1)) / kKernelNorm);
}

template <class TCallback>
void RunPassInParallel(int actual_threads, int stripes_per_worker, int extra_stripes, TCallback callback) {
  std::vector<std::thread> workers;
  workers.reserve(actual_threads);

  int stripe_begin = 0;
  for (int worker_id = 0; worker_id < actual_threads; ++worker_id) {
    const int stripes_this_worker = stripes_per_worker + (worker_id < extra_stripes ? 1 : 0);
    const int stripe_end = stripe_begin + stripes_this_worker;
    workers.emplace_back([callback, stripe_begin, stripe_end]() { callback(stripe_begin, stripe_end); });
    stripe_begin = stripe_end;
  }

  for (auto &worker : workers) {
    worker.join();
  }
}
}  // namespace

PikhotskiyRVerticalGaussFilterALL::PikhotskiyRVerticalGaussFilterALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool PikhotskiyRVerticalGaussFilterALL::ValidationImpl() {
  const auto &in = GetInput();

  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  const auto expected_size = static_cast<std::size_t>(in.width) * static_cast<std::size_t>(in.height);
  return in.data.size() == expected_size;
}

bool PikhotskiyRVerticalGaussFilterALL::PreProcessingImpl() {
  const auto &in = GetInput();
  width_ = in.width;
  height_ = in.height;

  const int num_threads = std::max(1, ppc::util::GetNumThreads());
  stripe_width_ = std::max(1, width_ / num_threads);

  source_ = in.data;
  vertical_buffer_.assign(source_.size(), 0);
  result_buffer_.assign(source_.size(), 0);
  return true;
}

bool PikhotskiyRVerticalGaussFilterALL::RunImpl() {
  const auto expected_size = static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_);
  if (width_ <= 0 || height_ <= 0 || source_.size() != expected_size || vertical_buffer_.size() != expected_size ||
      result_buffer_.size() != expected_size) {
    return false;
  }

  const int requested_threads = std::max(1, ppc::util::GetNumThreads());
  const int stripe_count = (width_ + stripe_width_ - 1) / stripe_width_;
  const int actual_threads = std::max(1, std::min(requested_threads, stripe_count));
  const int stripes_per_worker = stripe_count / actual_threads;
  const int extra_stripes = stripe_count % actual_threads;
  RunPassInParallel(actual_threads, stripes_per_worker, extra_stripes, [this](int stripe_begin, int stripe_end) {
    for (int stripe_index = stripe_begin; stripe_index < stripe_end; ++stripe_index) {
      const int x_begin = stripe_index * stripe_width_;
      const int x_end = std::min(width_, x_begin + stripe_width_);
      RunVerticalPassForStripe(x_begin, x_end);
    }
  });

  RunPassInParallel(actual_threads, stripes_per_worker, extra_stripes, [this](int stripe_begin, int stripe_end) {
    for (int stripe_index = stripe_begin; stripe_index < stripe_end; ++stripe_index) {
      const int x_begin = stripe_index * stripe_width_;
      const int x_end = std::min(width_, x_begin + stripe_width_);
      RunHorizontalPassForStripe(x_begin, x_end);
    }
  });

  return true;
}

bool PikhotskiyRVerticalGaussFilterALL::PostProcessingImpl() {
  GetOutput().width = width_;
  GetOutput().height = height_;
  GetOutput().data = result_buffer_;
  return true;
}

void PikhotskiyRVerticalGaussFilterALL::RunVerticalPassForStripe(int x_begin, int x_end) {
  for (int row = 0; row < height_; ++row) {
    const int row_top = ClampIndex(row - 1, height_);
    const int row_bottom = ClampIndex(row + 1, height_);

    for (int col = x_begin; col < x_end; ++col) {
      const std::size_t center = ToLinearIndex(col, row, width_);
      const std::size_t top = ToLinearIndex(col, row_top, width_);
      const std::size_t bottom = ToLinearIndex(col, row_bottom, width_);
      vertical_buffer_[center] =
          static_cast<int>(source_[top]) + (2 * static_cast<int>(source_[center])) + static_cast<int>(source_[bottom]);
    }
  }
}

void PikhotskiyRVerticalGaussFilterALL::RunHorizontalPassForStripe(int x_begin, int x_end) {
  for (int row = 0; row < height_; ++row) {
    for (int col = x_begin; col < x_end; ++col) {
      const int col_left = ClampIndex(col - 1, width_);
      const int col_right = ClampIndex(col + 1, width_);
      const std::size_t center = ToLinearIndex(col, row, width_);
      const std::size_t left = ToLinearIndex(col_left, row, width_);
      const std::size_t right = ToLinearIndex(col_right, row, width_);
      const int weighted_sum = vertical_buffer_[left] + (2 * vertical_buffer_[center]) + vertical_buffer_[right];
      result_buffer_[center] = NormalizeAndRoundUp(weighted_sum);
    }
  }
}

}  // namespace pikhotskiy_r_vertical_gauss_filter
