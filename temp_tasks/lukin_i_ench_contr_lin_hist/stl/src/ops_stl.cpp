#include "lukin_i_ench_contr_lin_hist/stl/include/ops_stl.hpp"

#include <algorithm>
#include <thread>
#include <vector>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"

namespace lukin_i_ench_contr_lin_hist {

LukinITestTaskSTL::LukinITestTaskSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType(GetInput().size());
}

bool LukinITestTaskSTL::ValidationImpl() {
  return !(GetInput().empty());
}

bool LukinITestTaskSTL::PreProcessingImpl() {
  return true;
}

bool LukinITestTaskSTL::RunImpl() {
  const InType &input = GetInput();

  size_ = static_cast<int>(input.size());

  thread_count_ = static_cast<int>(std::thread::hardware_concurrency());
  chunk_size_ = static_cast<int>(input.size()) / thread_count_;

  std::vector<unsigned char> loc_mins(thread_count_, 255);
  std::vector<unsigned char> loc_maxs(thread_count_, 0);

  GetLocMinMax(loc_mins, loc_maxs);

  const int min = *std::ranges::min_element(loc_mins);
  const int max = *std::ranges::max_element(loc_maxs);

  if (max == min) {
    GetOutput() = GetInput();
    return true;
  }

  Process(min, max);

  return true;
}

bool LukinITestTaskSTL::PostProcessingImpl() {
  return true;
}

void LukinITestTaskSTL::GetLocMinMax(std::vector<unsigned char> &loc_mins, std::vector<unsigned char> &loc_maxs) {
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

void LukinITestTaskSTL::Process(int min, int max) {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  const float scale = 255.0F / static_cast<float>(max - min);

  auto process = [&](const int start, const int end) {
    for (int i = start; i < end; i++) {
      output[i] = static_cast<unsigned char>(static_cast<float>(input[i] - min) * scale);
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

}  // namespace lukin_i_ench_contr_lin_hist
