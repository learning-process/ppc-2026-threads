#include <functional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace vinyaikina_e_multidimensional_integrals_simpson_method {

using InType =
    std::tuple<double, std::vector<std::pair<double, double>>, std::function<double(const std::vector<double> &)>>;
using OutType = double;

using TestType = std::tuple<std::string, InType, double>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
