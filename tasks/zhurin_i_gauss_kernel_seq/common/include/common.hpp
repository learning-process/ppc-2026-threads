#pragma once

#include <tuple>
#include <vector>

namespace zhurin_i_gauss_kernel_seq {

// Входные данные: (ширина, высота, количество вертикальных полос, изображение)
using InType = std::tuple<int, int, int, std::vector<std::vector<int>>>;

// Выходные данные: отфильтрованное изображение той же размерности
using OutType = std::vector<std::vector<int>>;

}  // namespace zhurin_i_gauss_kernel_seq