#pragma once

#include <functional>
#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace vlasova_a_simpson_method_seq {

struct SimpsonTask {
  std::function<double(const std::vector<double> &)> func;
  std::vector<double> a;
  std::vector<double> b;
  std::vector<int> n;

  SimpsonTask() = default;
  SimpsonTask(std::function<double(const std::vector<double> &)> f, const std::vector<double> &lower,
              const std::vector<double> &upper, const std::vector<int> &steps)
      : func(f), a(lower), b(upper), n(steps) {}
};

using InType = SimpsonTask;
using OutType = double;
using TestType = std::tuple<std::vector<double>, std::vector<double>, std::vector<int>, double, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace vlasova_a_simpson_method_seq
