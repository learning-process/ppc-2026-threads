#include <omp.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "savva_d_monte_carlo/common/include/common.hpp"
#include "savva_d_monte_carlo/omp/include/ops_omp.hpp"

namespace savva_d_monte_carlo {

SavvaDMonteCarloOMP::SavvaDMonteCarloOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SavvaDMonteCarloOMP::ValidationImpl() {
  const auto &input = GetInput();

  // Проверка количества точек
  if (input.count_points == 0) {
    return false;
  }

  // Проверка наличия функции
  if (!input.f) {
    return false;
  }

  // Проверка размерности
  if (input.Dimension() == 0) {
    return false;
  }

  // Проверка корректности границ
  for (size_t i = 0; i < input.Dimension(); ++i) {
    if (input.lower_bounds[i] > input.upper_bounds[i]) {
      return false;
    }
  }

  return true;
}

bool SavvaDMonteCarloOMP::PreProcessingImpl() {
  return true;
}

bool SavvaDMonteCarloOMP::RunImpl() {
  const auto &input = GetInput();
  auto &result = GetOutput();

  const size_t dim = input.Dimension();
  const double vol = input.Volume();
  // OpenMP исторически лучше работает со знаковыми типами счетчиков циклов
  const int64_t n = static_cast<int64_t>(input.count_points);
  const auto &func = input.f;

  std::vector<std::uniform_real_distribution<double>> distributions(dim);
  for (size_t i = 0; i < dim; ++i) {
    distributions[i] = std::uniform_real_distribution<double>(input.lower_bounds[i], input.upper_bounds[i]);
  }

  double sum = 0.0;

  // Вычисляем количество полных блоков по 4 элемента и остаток
  const int64_t n_blocks = n / 4;
  const int64_t tail = n % 4;

// Открываем параллельную секцию. Суммирование собираем через редукцию
#pragma omp parallel reduction(+ : sum)
  {
    // Каждый поток инициализирует свой собственный генератор.
    // XOR с номером потока гарантирует уникальность последовательностей.
    std::minstd_rand generator(std::random_device{}() ^ omp_get_thread_num());

    // Выделяем память под координаты один раз для каждого потока!
    // Это избавляет от накладных расходов на аллокацию внутри цикла.
    std::vector<double> p1(dim);
    std::vector<double> p2(dim);
    std::vector<double> p3(dim);
    std::vector<double> p4(dim);

// Параллельный цикл по блокам
#pragma omp for schedule(static)
    for (int64_t i = 0; i < n_blocks; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        p1[j] = distributions[j](generator);
        p2[j] = distributions[j](generator);
        p3[j] = distributions[j](generator);
        p4[j] = distributions[j](generator);
      }
      sum += func(p1) + func(p2) + func(p3) + func(p4);
    }
  }

  // Обрабатываем хвост последовательно (максимум 3 итерации, параллелить нет смысла)
  if (tail > 0) {
    std::minstd_rand generator(std::random_device{}());
    std::vector<double> p_tail(dim);
    for (int64_t i = 0; i < tail; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        p_tail[j] = distributions[j](generator);
      }
      sum += func(p_tail);
    }
  }

  // Вычисляем интеграл: среднее значение функции умноженное на объем области
  double mean = sum / static_cast<double>(n);
  result = mean * vol;

  return true;
}

bool SavvaDMonteCarloOMP::PostProcessingImpl() {
  return true;
}

}  // namespace savva_d_monte_carlo
