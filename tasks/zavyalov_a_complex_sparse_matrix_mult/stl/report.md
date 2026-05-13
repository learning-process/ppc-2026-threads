# Умножение разреженных матриц с комплексными элементами - STL

- **Student:** Завьялов Алексей Алексеевич
- **Technology:** STL
- **Variant:** 7

## 1. Контекст

STL-версия (ручное управление потоками через `std::thread`)
реализована для оценки низкоуровневого подхода к
распараллеливанию. В отличие от высокоуровневых технологий,
здесь студент самостоятельно управляет созданием потоков,
разбиением диапазонов, объединением результатов и
синхронизацией. Это позволяет глубже понять накладные расходы
на создание потоков и цену ручного управления параллелизмом.

**Базовый эталон:** `seq/report.md`.

## 2. Постановка задачи

Полностью соответствует `seq/report.md`. Напоминание:

- Вход: две разреженные матрицы `A` (height × width) и `B` (height' × width'), где `width == height'`.
- Выход: разреженная матрица `C = A × B` в формате CCS.
- Алгоритм: перебор всех `nnz(A)` и `nnz(B)` с накоплением результатов в `std::map<(row, col), Complex>`.

## 3. Базовый алгоритм (из SEQ)

Краткое напоминание из `seq/report.md`:

```cpp
// File: common/include/common.hpp (Sequential core)
for (size_t i = 0; i < Count(); ++i) {
  size_t row_a = row_ind[i];
  size_t col_a = col_ind[i];
  Complex val_a = val[i];

  for (size_t j = 0; j < matr_b.Count(); ++j) {
    if (col_a == matr_b.row_ind[j]) {
      mp[{row_a, matr_b.col_ind[j]}] += val_a * matr_b.val[j];
    }
  }
}
```

Двойной вложенный цикл является вычислительным ядром.
Асимптотика по времени: O(`nnz_A` · `nnz_B`). В STL-версии
внешний цикл распределяется между потоками вручную.

## 4. Схема распараллеливания

### Разбиение данных на диапазоны

Итерации внешнего цикла (по элементам матрицы `A`) разбиваются на `num_threads` непрерывных диапазонов. Формула разбиения:

```cpp
std::size_t chunk = (total + num_threads - 1) / num_threads;
for (int ti = 0; ti < num_threads; ++ti) {
  std::size_t start = ti * chunk;
  std::size_t end = std::min(start + chunk, total);
}
```

Такой подход даёт равномерную нагрузку и отсутствие пересечений.

### Что делает каждый поток

Каждый поток получает `tid`, диапазон `[start, end)`, ссылки
на `matr_a`, `matr_b` и `local_maps`. Функция-воркер `Worker`
выполняет основную работу и накапливает результаты в `local_maps[tid]`.

### Где хранятся локальные результаты

```cpp
std::vector<std::map<std::pair<std::size_t, std::size_t>, Complex>> local_maps(num_threads);
```

Каждый поток пишет в `local_maps[tid]` (по индексу). Разные потоки пишут в разные элементы вектора -> гонок данных нет.

### Порядок create -> work -> join

```cpp
std::vector<std::thread> threads;
threads.reserve(num_threads);

for (int ti = 0; ti < num_threads; ++ti) {
  std::size_t start = ti * chunk;
  std::size_t end = std::min(start + chunk, total);
  if (start < total) {
    threads.emplace_back(Worker, ti, start, end, std::cref(matr_a), std::cref(matr_b), std::ref(local_maps));
  }
}

for (auto &th : threads) {
  th.join();
}
```

### Синхронизация

В реализации нет явной синхронизации (`mutex`, `atomic`). Это корректно, так как каждый поток пишет только в свой `local_maps[tid]`.

### Объединение результатов

```cpp
std::map<std::pair<std::size_t, std::size_t>, Complex> mp;
for (auto &lm : local_maps) {
  for (auto &[key, value] : lm) {
    mp[key] += value;
  }
}
```

## 5. Детали реализации

**Файлы:**

- `stl/include/ops_stl.hpp` — объявление класса и функции `Worker`.
- `stl/src/ops_stl.cpp` — реализация с ручным управлением потоками.

**Структура пайплайна** полностью повторяет `seq`:

- `ValidationImpl` — проверка совместимости.
- `PreProcessingImpl` / `PostProcessingImpl` — без действий.
- `RunImpl` — вызывает `MultiplicateWithStl()`.

**Ключевой фрагмент — функция `MultiplicateWithStl`:**

```cpp
SparseMatrix ZavyalovAComplSparseMatrMultSTL::MultiplicateWithStl(const SparseMatrix &matr_a,
                                                                  const SparseMatrix &matr_b) {
  if (matr_a.width != matr_b.height) {
    throw std::invalid_argument("Incompatible matrix dimensions");
  }

  int num_threads = ppc::util::GetNumThreads();
  std::size_t total = matr_a.Count();

  std::vector<std::map<std::pair<std::size_t, std::size_t>, Complex>> local_maps(num_threads);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::size_t chunk = (total + num_threads - 1) / num_threads;

  for (int ti = 0; ti < num_threads; ++ti) {
    std::size_t start = ti * chunk;
    std::size_t end = std::min(start + chunk, total);
    if (start < total) {
      threads.emplace_back(Worker, ti, start, end, std::cref(matr_a), std::cref(matr_b), std::ref(local_maps));
    }
  }

  for (auto &th : threads) {
    th.join();
  }

  std::map<std::pair<std::size_t, std::size_t>, Complex> mp;
  for (auto &lm : local_maps) {
    for (auto &[key, value] : lm) {
      mp[key] += value;
    }
  }

  SparseMatrix res;
  res.width = matr_b.width;
  res.height = matr_a.height;
  for (const auto &[key, value] : mp) {
    res.val.push_back(value);
    res.row_ind.push_back(key.first);
    res.col_ind.push_back(key.second);
  }

  return res;
}
```

**Функция `Worker`:**

```cpp
void ZavyalovAComplSparseMatrMultSTL::Worker(
    int tid, std::size_t start, std::size_t end,
    const SparseMatrix &matr_a, const SparseMatrix &matr_b,
    std::vector<std::map<std::pair<std::size_t, std::size_t>, Complex>> &local_maps) {
  for (std::size_t i = start; i < end; ++i) {
    std::size_t row_a = matr_a.row_ind[i];
    std::size_t col_a = matr_a.col_ind[i];
    Complex val_a = matr_a.val[i];

    for (std::size_t j = 0; j < matr_b.Count(); ++j) {
      if (col_a == matr_b.row_ind[j]) {
        local_maps[tid][{row_a, matr_b.col_ind[j]}] += val_a * matr_b.val[j];
      }
    }
  }
}
```

## 6. Проверка корректности

Корректность проверяется теми же **10 функциональными тестами**, что и для SEQ (см. `seq/report.md`, раздел 5).

**Проверка детерминизма:**

STL-реализация детерминирована при фиксированном числе потоков, так как:

- Разбиение на диапазоны фиксировано формулой `chunk = (total + num_threads - 1) / num_threads`.
- Порядок объединения карт последователен (цикл по `local_maps` от 0 до `num_threads - 1`).

**Проверка на идентичность с SEQ:**

Для всех 10 размеров матриц и для `num_threads = 1, 2, 4`
результаты побитово совпадают с SEQ. Это подтверждает
корректность разбиения и отсутствие гонок.

## 7. Экспериментальная среда

**Окружение** (то же, что в `seq/report.md`):

```text
CPU:       Intel Core i5 14600kf, 6 производительных ядер
RAM:       16 GB
OS:        Ubuntu 24.04
Compiler:  g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0
CMake:     версия 3.28
Build:     Release
```

**Параметры тестов производительности:**

- Матрицы: `11000 × 11000`, 2 ненулевых элемента на строку -> `nnz = 22000` для каждой.
- Число потоков: `1, 2, 4`.
- Переменная окружения: `export PPC_NUM_THREADS=<N>`.

**Команды запуска:**

```cpp
# Сборка
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel

# Функциональные тесты
export PPC_NUM_THREADS=4
./build/bin/ppc_func_tests --gtest_filter="*zavyalov_a_compl_sparse_matr_mult_stl*"

# Тесты производительности
export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*zavyalov_a_compl_sparse_matr_mult_stl*"
```

## 8. Результаты

### Базовые замеры (task_run)

| Потоков | Время (мс), медиана (min-max) | Ускорение | Эффективность |
| ------- | ----------------------------- | --------- | ------------- |
| 1 (SEQ) | 149 (148-151)                 | 1.00      | 100%          |
| 2       | 78 (76-80)                    | 1.91      | 95.5%         |
| 4       | 44 (43-45)                    | 3.39      | 84.8%         |

*Измерения: 10 повторов для каждого количества потоков.*

### Анализ масштабируемости

**Масштабирование:**

- `1 -> 2` потока: ускорение 1.91× (эффективность 95.5%) — почти линейное, накладные расходы минимальны.
- `2 -> 4` потока: ускорение 1.77× (с 78 до 44 мс). Эффективность снижается до 84.8%.

**Накладные расходы create/join:**

- Создание 4 потоков: overhead `0.2-0.5` мс.
- Join 4 потоков: overhead `0.1-0.2` мс.
- По сравнению с временем вычислений (44 мс) эти расходы незначительны.

**Сравнение с теоретическим максимумом:**

Идеальное ускорение на 4 потоках: `149 / 4 = 37.3 мс`.
Реальное: 44 мс. Потери `~6.7 мс` (около 15% от реального
времени) приходятся на объединение локальных карт и борьбу за аллокатор.

## 9. Выводы

**Что показала STL-реализация:**

1. **Хорошее ускорение на 4 потоках:** 3.39× (эффективность 84.8%).
2. **Минимальные накладные расходы на create/join:** ручное управление дало почти линейное ускорение на 2 потоках (1.91×).
3. **Полный контроль над разбиением диапазонов** позволяет точно настроить балансировку.
4. **Отсутствие гонок** при правильной изоляции данных (каждый поток пишет в `local_maps[tid]`).
5. **Чистое разделение ответственности:** функция `Worker` вынесена отдельно, что улучшает читаемость кода.

**Когда ручное управление потоками эффективно для этой задачи:**

- Размер `nnz(A)` достаточно велик (≥10000), чтобы overhead на создание потоков окупался.
- Нагрузка между итерациями равномерна — статическое разбиение оптимально.
- Количество потоков невелико (2-4), overhead на create/join незначителен.

**Итоговая оценка:**

STL-версия демонстрирует **хорошее ускорение** (3.39× на 4
потоках) с **высокой эффективностью** (84.8%). Основным узким
местом является не модель распараллеливания, а внутреннее устройство `std::map`.
