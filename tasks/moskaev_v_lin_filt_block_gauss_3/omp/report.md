# Линейная фильтрация блоков – фильтр Гаусса 3×3 — OMP

- Student: Москаев Владимир Александрович
- Technology: OMP
- Variant: 26

## 1. Контекст

OpenMP-реализация блочной фильтрации изображения ядром Гаусса 3×3.
Распараллелены три уровня: фильтрация внутри блока (collapse(3)),
копирование с отступами (collapse(2)) и внешний цикл по блокам.
Цель — оценить эффективность OpenMP на задаче с регулярной вычислительной нагрузкой.

## 2. Постановка задачи

Совпадает с последовательной версией (см. seq/report.md).
Вход — изображение, выход — отфильтрованное изображение.
Корректность проверяется сравнением с SEQ.

## 3. Базовый алгоритм

Тот же блочный алгоритм, что в SEQ: разбиение на блоки 64×64,
копирование с отступами, свёртка с ядром Гаусса, запись результата.

## 4. Схема распараллеливания

**Какая область кода параллелится:**

- Фильтрация внутри блока (тройной вложенный цикл)
- Копирование блока с отступами (двойной цикл)
- Внешний цикл по блокам (двойной цикл по y и x)

**Какие переменные shared/private:**

- shared: input_block, output_block, inner_width, inner_height, channels,
  block_width, kGaussianKernel, image_data, width, height, block_size
- private: row, col, channel, ky, kx, sum, idx, out_idx (неявно через collapse)

**Где нужен reduction / atomic / critical:**
Не требуется. Каждый поток пишет в свои независимые выходные данные
(разные блоки или разные пиксели внутри блока). Гонок нет.

**Какой schedule выбран и почему:**

- Для фильтрации и копирования: schedule(static) — минимальный оверхед,
  нагрузка равномерна.
- Для внешнего цикла по блокам: schedule(dynamic) — блоки у границ могут быть меньше,
  dynamic лучше балансирует.

**Неявный барьер:**
В конце каждой параллельной области `#pragma omp parallel for` существует неявный барьер.
Это означает, что все потоки синхронизируются перед тем,
как выполнение продолжится на главном потоке.
В данной реализации барьер не влияет на производительность,
так как после фильтрации каждого блока нет зависимостей по данным между итерациями.

## 5. Детали реализации

**Файлы:** omp/include/ops_omp.hpp, omp/src/ops_omp.cpp

**Какие участки были изменены относительно SEQ:**

- В ApplyGaussianFilterToBlock() добавлена директива:

```cpp
  #pragma omp parallel for collapse(3) schedule(static) default(none) shared(...)
```

- В CopyBlockWithPadding() добавлена:

```cpp
  #pragma omp parallel for collapse(2) schedule(static) default(none) shared(...)
  ```

- В CopyProcessedBlockToOutput() добавлена:

```cpp
  #pragma omp parallel for collapse(2) schedule(static) default(none) shared(...)
  ```

- В RunImpl() внешний цикл по блокам:

```cpp
  #pragma omp parallel for collapse(2) schedule(dynamic) default(none) shared(...)
  ```

**Фрагмент кода фильтрации:**

```cpp
#pragma omp parallel for collapse(3) schedule(static) default(none) \
    shared(input_block, output_block, inner_width, inner_height, channels, block_width, kGaussianKernel)
for (int row = 0; row < inner_height; ++row) {
    for (int col = 0; col < inner_width; ++col) {
        for (int channel = 0; channel < channels; ++channel) {
            float sum = 0.0f;
            for (int ky = -1; ky <= 1; ++ky)
                for (int kx = -1; kx <= 1; ++kx)
                    sum += input_block[...] * kGaussianKernel[...];
            output_block[...] = static_cast<uint8_t>(std::round(sum));
        }
    }
}
```

**Какие риски гонок были устранены:**
Гонок нет, так как:

- При фильтрации каждый поток работает с разными пикселями (итерации цикла независимы).
- При копировании блоков каждый поток обрабатывает разные блоки.

## 6. Проверка корректности

Сравнение с SEQ. Функциональные тесты (из tests/functional/main.cpp):

- Тест 1: 2×2 серое, вход [100,150,200,250] → выход [138,163,188,213]
- Тест 2: 3×3 серое, вход [1,2,3,4,5,6,7,8,9] → выход [2,3,4,4,5,6,7,7,8]
- Тест 3: 2×2 RGB, вход 12 чисел → выход 12 чисел
- Тест 4: 1×1 серое, вход [255] → выход [255]

Все тесты пройдены для 1, 2, 4, 8 потоков.

## 7. Экспериментальная среда

- CPU: Intel Core i5-11400H
- RAM: 16 ГБ
- OS: Windows 10
- Компилятор: MSVC 2022 с поддержкой OpenMP
- Размер задачи: 2048×2048×3 пикселя

**Переменные окружения:** OMP_NUM_THREADS

**Команды запуска:**

```bash
$env:OMP_NUM_THREADS=1; .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_moskaev_v_lin_filt_block_gauss_3_omp_enabled"
$env:OMP_NUM_THREADS=2; .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_moskaev_v_lin_filt_block_gauss_3_omp_enabled"
$env:OMP_NUM_THREADS=4; .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_moskaev_v_lin_filt_block_gauss_3_omp_enabled"
$env:OMP_NUM_THREADS=8; .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_moskaev_v_lin_filt_block_gauss_3_omp_enabled"
```

## 8. Результаты

| Потоков | Время, с | Speedup | Efficiency |
|---------|----------|---------|------------|
| 1       | 0.645    | 0.33    | 33%        |
| 2       | 0.293    | 0.73    | 37%        |
| 4       | 0.161    | 1.33    | 33%        |
| 8       | 0.115    | 1.86    | 23%        |

**Комментарий о масштабируемости и узких местах:**

- При 1 потоке OpenMP в 3 раза медленнее SEQ из-за оверхеда на создание параллельных областей.
- Ускорение растёт до 1.86× на 8 потоках, но далеко от идеального (8×).
- Узкое место: последовательная природа внешнего цикла (каждый шаг зависит от предыдущего) ограничивает максимальное ускорение.
- Эффективность падает с ростом числа потоков из-за увеличения доли оверхеда.

## 9. Выводы

OpenMP даёт ускорение до 1.86× на 8 потоках.
При малом числе потоков (1-2) оверхед перевешивает выигрыш.
Рекомендуется использовать не менее 4 потоков.
