#include "zaharov_g_linear_contrast_stretch/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <omp.h>

#include "util/include/util.hpp"
#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"

namespace zaharov_g_linear_contrast_stretch {

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
  int min_el = 255;
  int max_el = 0;

  for (std::int64_t i = 0; i < size; ++i) {
    const int value = static_cast<int>(input[static_cast<size_t>(i)]);
    min_el = std::min(min_el, value);
    max_el = std::max(max_el, value);
  }

  if (max_el > min_el) {
    const int denom = max_el - min_el;

    for (std::int64_t i = 0; i < size; ++i) {
      const int value = (static_cast<int>(input[static_cast<size_t>(i)]) - min_el) * 255 / denom;
      output[static_cast<size_t>(i)] = static_cast<uint8_t>(std::clamp(value, 0, 255));
    }
  } else {
    for (std::int64_t i = 0; i < size; ++i) {
      output[static_cast<size_t>(i)] = input[static_cast<size_t>(i)];
    }
  }

  return true;
}

bool ZaharovGLinContrStrOMP::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace zaharov_g_linear_contrast_stretch
