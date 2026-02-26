#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kondrashova_v_marking_components_seq {

struct ImageData {
  std::vector<uint8_t> data;
  int width{};
  int height{};
};

struct Result {
  int count{};
  std::vector<std::vector<int>> labels;
};

using InType = ImageData;
using OutType = Result;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kondrashova_v_marking_components_seq
