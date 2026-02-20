#include "lukin_i_ench_contr_lin_hist/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"
#include "util/include/util.hpp"

namespace lukin_i_ench_contr_lin_hist {

LukinITestTaskSEQ::LukinITestTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType(GetInput().size());
}

bool LukinITestTaskSEQ::ValidationImpl() {
  return !(GetInput().empty());
}

bool LukinITestTaskSEQ::PreProcessingImpl() {
  return true;
}

bool LukinITestTaskSEQ::RunImpl() {
  unsigned char min = 255;
  unsigned char max = 0;

  for (const auto &elem : GetInput())  // Поиск максимума и минимума
  {
    if (elem > max) {
      max = elem;
    }
    if (elem < min) {
      min = elem;
    }
  }

  if (max == min)  // Однотонное изображение
  {
    GetOutput() = GetInput();
    return true;
  }

  float scale = 255.0f / (max - min);
  for (int i = 0; i < static_cast<int>(GetInput().size()); i++) {  // Линейное растяжение
    GetOutput()[i] = static_cast<unsigned char>((GetInput()[i] - min) * scale);
  }

  return true;
}

bool LukinITestTaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace lukin_i_ench_contr_lin_hist
