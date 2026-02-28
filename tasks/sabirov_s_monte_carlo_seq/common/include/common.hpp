#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace sabirov_s_monte_carlo_seq {

enum class FuncType : uint8_t {
  kLinear = 0,
  kSumCubes = 1,
  kCosProduct = 2,
  kExpNeg = 3,
  kMixedPoly = 4,
  kSinSum = 5,
  kSqrtSum = 6,
  kQuarticSum = 7,
};

struct MCInput {
  std::vector<double> lower;
  std::vector<double> upper;
  int num_samples;
  FuncType func_type;
};

using InType = MCInput;
using OutType = double;
using TestType = std::tuple<InType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sabirov_s_monte_carlo_seq
