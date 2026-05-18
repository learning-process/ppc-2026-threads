#include "pylaeva_s_inc_contrast_img_by_lsh/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

PylaevaSIncContrastImgByLshSTL::PylaevaSIncContrastImgByLshSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType(GetInput().size());
}

bool PylaevaSIncContrastImgByLshSTL::ValidationImpl() {
  return !(GetInput().empty());
}

bool PylaevaSIncContrastImgByLshSTL::PreProcessingImpl() {
  return true;
}

bool PylaevaSIncContrastImgByLshSTL::RunImpl() {
  const InType &input = GetInput();

  size_ = static_cast<int>(input.size());

  thread_count_ = static_cast<int>(std::thread::hardware_concurrency());
  chunk_size_ = static_cast<int>(input.size()) / thread_count_;

  std::vector<unsigned char> loc_mins(thread_count_, 255);
  std::vector<unsigned char> loc_maxs(thread_count_, 0);

  FindLocalMinMax(loc_mins, loc_maxs);

  const int min = *std::ranges::min_element(loc_mins);
  const int max = *std::ranges::max_element(loc_maxs);

  if (max == min) {
    GetOutput() = GetInput();
    return true;
  }

  ApplyLinearStretching(min, max);

  return true;
}

bool PylaevaSIncContrastImgByLshSTL::PostProcessingImpl() {
  return true;
}

void PylaevaSIncContrastImgByLshSTL::FindLocalMinMax(std::vector<unsigned char> &loc_mins,
                                                     std::vector<unsigned char> &loc_maxs) {
  const InType &input = GetInput();

  std::vector<std::thread> thread_pool;
  thread_pool.reserve(thread_count_);

  auto reduction = [&](const int idx, const int start, const int end) {
    unsigned char loc_min = loc_mins[idx];
    unsigned char loc_max = loc_maxs[idx];

    for (int i = start; i < end; i++) {
      loc_min = (input[i] < loc_min) ? input[i] : loc_min;
      loc_max = (input[i] > loc_max) ? input[i] : loc_max;
    }

    loc_mins[idx] = loc_min;
    loc_maxs[idx] = loc_max;
  };

  for (int i = 0; i < thread_count_; i++) {
    const int start = chunk_size_ * i;
    const int end = (i == (thread_count_ - 1)) ? size_ : start + chunk_size_;
    thread_pool.emplace_back(reduction, i, start, end);
  }
  for (auto &thread : thread_pool) {
    thread.join();
  }
}

void PylaevaSIncContrastImgByLshSTL::ApplyLinearStretching(int min, int max) {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  const double scale = 255.0 / static_cast<double>(max - min);
  const double offset = static_cast<double>(-min) * scale;

  auto process = [&](const int start, const int end) {
    for (int i = start; i < end; i++) {
      // Предварительно вычисленное выражение: input[i] * scale + offset
      double new_value = (static_cast<double>(input[i]) * scale) + offset;
      // Округление к ближайшему целому
      int rounded_value = static_cast<int>(std::lround(new_value));

      rounded_value = std::clamp(rounded_value, 0, 255);

      output[i] = static_cast<unsigned char>(rounded_value);
    }
  };

  std::vector<std::thread> thread_pool;
  thread_pool.reserve(thread_count_);

  for (int i = 0; i < thread_count_; i++) {
    const int start = chunk_size_ * i;
    const int end = (i == (thread_count_ - 1)) ? size_ : start + chunk_size_;
    thread_pool.emplace_back(process, start, end);
  }
  for (auto &thread : thread_pool) {
    thread.join();
  }
}

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
