#include "zaharov_g_linear_contrast_stretch/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"
#include "util/include/util.hpp"

namespace zaharov_g_linear_contrast_stretch {

ZaharovGLinContrStrSEQ::ZaharovGLinContrStrSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZaharovGLinContrStrSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool ZaharovGLinContrStrSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool ZaharovGLinContrStrSEQ::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  auto [min_it, max_it] = std::minmax_element(input.begin(), input.end());
  uint8_t min_el = *min_it, max_el = *max_it;

  if (max_el > min_el) {
    int denom = max_el - min_el;

    for (size_t i = 0; i < input.size(); ++i) {
      int value = (static_cast<int>(input[i]) - min_el) * 255 / denom;
      output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
    }
  } else {
    output.assign(input.begin(), input.end());
  }

  return true;
}

bool ZaharovGLinContrStrSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace zaharov_g_linear_contrast_stretch
