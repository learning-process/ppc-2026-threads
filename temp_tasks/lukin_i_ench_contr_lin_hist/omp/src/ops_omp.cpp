#include "lukin_i_ench_contr_lin_hist/omp/include/ops_omp.hpp"

#include <algorithm>
#include <vector>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"

namespace lukin_i_ench_contr_lin_hist {

LukinITestTaskOMP::LukinITestTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType(GetInput().size());
}

bool LukinITestTaskOMP::ValidationImpl() {
  return !(GetInput().empty());
}

bool LukinITestTaskOMP::PreProcessingImpl() {
  return true;
}

bool LukinITestTaskOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  unsigned char min = 255;
  unsigned char max = 0;

  const int size = static_cast<int>(input.size());

#pragma omp parallel for default(none) shared(input, size) reduction(min : min) \
    reduction(max : max)  //(max : max) - оператор, потом указание переменной
  for (int i = 0; i < size; i++) {
    min = std::min(min, input[i]);
    max = std::max(max, input[i]);
  }

  if (max == min)  // Однотонное изображение
  {
    output = input;
    return true;
  }

  float scale = 255.0F / static_cast<float>(max - min);

#pragma omp parallel for default(none) shared(input, output, min, size, scale)
  for (int i = 0; i < size; i++) {  // Линейное растяжение
    output[i] = static_cast<unsigned char>(static_cast<float>(input[i] - min) * scale);
  }

  return true;
}

bool LukinITestTaskOMP::PostProcessingImpl() {
  return true;
}

}  // namespace lukin_i_ench_contr_lin_hist
