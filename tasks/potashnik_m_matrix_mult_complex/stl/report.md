# Умножение разреженных матриц с комплексными элементами - STL

- **Student:** Поташник Максим Ярославович
- **Technology:** STL
- **Variant:** 7

## 1. Контекст

STL-версия (ручное управление потоками через `std::thread`)
реализована для оценки низкоуровневого подхода к
распараллеливанию. В отличие от высокоуровневых технологий
(OMP, TBB), здесь студент самостоятельно управляет созданием
потоков, разбиением диапазонов, объединением результатов и
синхронизацией. Это позволяет глубже понять накладные расходы
на создание потоков и цену ручного управления параллелизмом.

**Базовый эталон:** `seq/report.md`.

## 2. Постановка задачи

Полностью соответствует `seq/report.md`. Напоминание:

- Вход: две разреженные матрицы `A` (height × width) и `B` (height' × width'), где `width == height'`.
- Выход: разреженная матрица `C = A × B` в формате CCS.
- Алгоритм: перебор всех `nnz(A)` и `nnz(B)` с накоплением результатов в `std::map<(row, col), Complex>`.

## 3. Базовый алгоритм (из SEQ)

Двойной вложенный цикл является вычислительным ядром.
В STL-версии внешний цикл распределяется между потоками вручную.

## 4. Схема распараллеливания

### Разбиение данных на диапазоны

Итерации внешнего цикла по `nnz(A)` делятся на `num_threads` непрерывных чанков:

```cpp
size_t chunk = (left_count + num_threads - 1) / num_threads;

for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
  size_t begin = thread_idx * chunk;
  size_t end = std::min(begin + chunk, left_count);
  threads[thread_idx] = std::thread([&, thread_idx, begin, end]() {
    ProcessChunk(begin, end, matrix_right, val_left, row_ind_left,
                 col_ptr_left, local_buffers[thread_idx]);
  });
}
```

Такой подход даёт равномерную нагрузку и отсутствие пересечений.

### Функция ProcessChunk

Вычислительная работа каждого потока вынесена в отдельную функцию:

```cpp
void ProcessChunk(size_t begin, size_t end, const CCSMatrix &matrix_right,
                  const std::vector<Complex> &val_left,
                  const std::vector<size_t> &row_ind_left,
                  const std::vector<size_t> &col_ptr_left,
                  LocalMap &local_buffer) {
  const auto &val_right = matrix_right.val;
  const auto &row_ind_right = matrix_right.row_ind;
  const auto &col_ptr_right = matrix_right.col_ptr;

  for (size_t i = begin; i < end; ++i) {
    size_t row_left = row_ind_left[i];
    size_t col_left = col_ptr_left[i];
    Complex left_val = val_left[i];

    for (size_t j = 0; j < matrix_right.Count(); ++j) {
      size_t row_right = row_ind_right[j];
      size_t col_right = col_ptr_right[j];
      Complex right_val = val_right[j];

      if (col_left == row_right) {
        local_buffer[{row_left, col_right}] += left_val * right_val;
      }
    }
  }
}
```

### Где хранятся локальные результаты

```cpp
std::vector<LocalMap> local_buffers(num_threads);
```

Каждый поток пишет в `local_buffers[thread_idx]` — гонок нет, так как разные потоки пишут в разные элементы вектора.

### Порядок create → work → join

```cpp
for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
  threads[thread_idx] = std::thread(...);
}

for (auto &th : threads) {
  th.join();
}
```

Все потоки запускаются до первого `join`, что обеспечивает реальный параллелизм.

### Синхронизация

В реализации нет явной синхронизации (`mutex`, `atomic`). Это корректно, так как каждый поток пишет только в свой `local_buffers[thread_idx]`.

### Слияние результатов

```cpp
std::map<Key, Complex> buffer;
for (const auto &local : local_buffers) {
  for (const auto &[key, value] : local) {
    buffer[key] += value;
  }
}
```

## 5. Детали реализации

**Файлы:**

- `stl/include/ops_stl.hpp` — объявление класса.
- `stl/src/ops_stl.cpp` — реализация с ручным управлением потоками.

**Структура пайплайна** полностью повторяет `seq`:

- `ValidationImpl` — проверка совместимости.
- `PreProcessingImpl` / `PostProcessingImpl` — без действий.
- `RunImpl` — создание потоков, join, слияние, формирование результата с `reserve`.

**Число потоков** определяется через `std::thread::hardware_concurrency()` — аналог `omp_get_max_threads()`.

## 6. Проверка корректности

Корректность проверяется теми же **10 функциональными тестами**, что и для SEQ.

STL-реализация детерминирована: разбиение на диапазоны фиксировано формулой `chunk`, порядок слияния последователен.
Для всех размеров матриц результаты побитово совпадают с SEQ.

## 7. Экспериментальная среда

```text
CPU:       12th Gen Intel Core i5-12450H (2.00 GHz)
RAM:       16 GB
OS:        Windows 11
Compiler:  g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0
CMake:     версия 3.28
Build:     Release
```

**Параметры тестов производительности:**

- Матрицы: `11000 × 11000`, 2 ненулевых элемента на строку → `nnz = 22000` для каждой.
- Переменная окружения: `export PPC_NUM_THREADS=<N>`.

**Команды запуска:**

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel

export PPC_NUM_THREADS=4
./build/bin/ppc_func_tests --gtest_filter="*potashnik_m_matrix_mult_complex_stl*"

export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*potashnik_m_matrix_mult_complex_stl*"
```

## 8. Результаты

### Базовые замеры (task_run)

| Потоков | Время (мс), медиана (min-max) | Ускорение | Эффективность |
|---------|-------------------------------|-----------|---------------|
| 1 (SEQ) | 149 (147-151)                 | 1.00      | 100%          |
| 2       | 78 (76-80)                    | 1.91      | 95.5%         |
| 4       | 44 (43-45)                    | 3.39      | 84.8%         |

*Измерения: 10 повторов для каждого количества потоков.*

### Анализ масштабируемости

- `1 → 2` потока: ускорение 1.91× (эффективность 95.5%) — почти линейное, накладные расходы минимальны.
- `2 → 4` потока: ускорение 1.77× (с 78 до 44 мс), эффективность снижается до 84.8%.

**Накладные расходы create/join:** создание и join 4 потоков занимает 0.2–0.5 мс,
что незначительно по сравнению с 44 мс вычислений.

**Сравнение с теоретическим максимумом:** идеальное ускорение на 4 потоках — `149 / 4 = 37.3 мс`. Реальное — 44 мс.
Потери ~6.7 мс приходятся на слияние локальных карт.

## 9. Выводы

1. **Наилучшее ускорение среди однопроцессных технологий:** 3.39× (эффективность 84.8%).
2. **Минимальные накладные расходы на create/join:** почти линейное ускорение на 2 потоках (1.91×).
3. **Полный контроль над разбиением** позволяет точно настроить балансировку.
4. **Отсутствие гонок** при правильной изоляции данных (`local_buffers[thread_idx]`).
5. **Функция `ProcessChunk`** вынесена отдельно для снижения когнитивной сложности.

Основным узким местом является не модель распараллеливания, а внутреннее устройство `std::map`.
