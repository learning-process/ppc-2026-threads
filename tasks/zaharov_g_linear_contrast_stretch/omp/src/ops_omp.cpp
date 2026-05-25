#include "zaharov_g_linear_contrast_stretch/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "util/include/util.hpp"
#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"

namespace zaharov_g_linear_contrast_stretch {

namespace {

struct MinMax {
  int min;
  int max;
};

MinMax FindMinMax(const InType &input) {
  const auto size = static_cast<std::int64_t>(input.size());
  const int thread_count = std::max(1, ppc::util::GetNumThreads());
  std::vector<MinMax> local_minmax(
      static_cast<std::size_t>(thread_count),
      MinMax{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()});

#pragma omp parallel default(none) shared(input, local_minmax, size) num_threads(thread_count)
  {
    MinMax current{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};
#pragma omp for nowait
    for (std::int64_t i = 0; i < size; ++i) {
      const int value = static_cast<int>(input[static_cast<std::size_t>(i)]);
      current.min = std::min(current.min, value);
      current.max = std::max(current.max, value);
    }
    local_minmax[static_cast<std::size_t>(omp_get_thread_num())] = current;
  }

  MinMax result{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};
  for (const auto &current : local_minmax) {
    result.min = std::min(result.min, current.min);
    result.max = std::max(result.max, current.max);
  }
  return result;
}

}  // namespace

ZaharovGLinContrStrOMP::ZaharovGLinContrStrOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZaharovGLinContrStrOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool ZaharovGLinContrStrOMP::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool ZaharovGLinContrStrOMP::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  const auto size = static_cast<std::int64_t>(input.size());
  const MinMax minmax = FindMinMax(input);

  if (minmax.max > minmax.min) {
    const int denom = minmax.max - minmax.min;

#pragma omp parallel for default(none) shared(input, output, size, minmax, denom) \
    num_threads(ppc::util::GetNumThreads())
    for (std::int64_t i = 0; i < size; ++i) {
      const int value = (static_cast<int>(input[static_cast<std::size_t>(i)]) - minmax.min) *
                        std::numeric_limits<uint8_t>::max() / denom;
      output[static_cast<std::size_t>(i)] = static_cast<uint8_t>(std::clamp(value, 0, 255));
    }
  } else {
#pragma omp parallel for default(none) shared(input, output, size) num_threads(ppc::util::GetNumThreads())
    for (std::int64_t i = 0; i < size; ++i) {
      output[static_cast<std::size_t>(i)] = input[static_cast<std::size_t>(i)];
    }
  }

  return true;
}

bool ZaharovGLinContrStrOMP::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace zaharov_g_linear_contrast_stretch
