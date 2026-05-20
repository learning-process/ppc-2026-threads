# Сортировка Хоара с простым слиянием — OMP

- **Студент:** Олесницкий Владимир Тарасович, 3823Б1ПР2
- **Технология:** OMP
- **Вариант:** 13

## 1. Контекст

OpenMP-версия расширяет базовую сортировку Хоара за счет параллельной обработки
независимых блоков и последующего простого слияния. Цель реализации — проверить,
окупаются ли параллельные циклы OpenMP на входе `N=100000`.

## 2. Постановка задачи

- **Входные данные:** непустой объект `std::vector<int>`.
- **Выходные данные:** отсортированный по неубыванию `std::vector<int>`.
- **Baseline:** последовательная версия с временем `T_seq = 0.0058254364 s`.

## 3. Базовый алгоритм

Последовательное ядро использует разбиение Хоара и quicksort по локальным
диапазонам. После локальной сортировки блоков по 64 элемента выполняется
попарное простое слияние до полного массива.

Средняя сложность quicksort-части составляет `O(n log n)`. Слияние добавляет
линейный проход на каждом уровне объединения блоков.

## 4. Схема распараллеливания

OpenMP используется в двух `parallel for`:

1. сортировка независимых блоков;
2. слияние независимых пар блоков.

В директивах используется `default(none)`, поэтому общие переменные перечислены
явно: `size`, `data`, `merged_data`, `merge_width`. Индексы циклов и временные
векторы являются локальными для тела цикла. Reduction не требуется, так как
общей скалярной агрегации нет.

Фрагмент параллельной части:

```cpp
#pragma omp parallel for default(none) shared(size, data)
  for (std::size_t block_start = 0; block_start < size;
       block_start += kBlockSize) {
    const std::size_t block_end = std::min(block_start + kBlockSize, size);
    if ((block_end - block_start) > 1) {
      HoareQuickSort(data, static_cast<int>(block_start),
                     static_cast<int>(block_end - 1));
    }
  }

  for (std::size_t merge_width = kBlockSize; merge_width < size;
       merge_width *= 2) {
    std::vector<int> merged_data(size);

#pragma omp parallel for default(none) \
    shared(merge_width, size, merged_data, data)
    for (std::size_t left = 0; left < size; left += (2 * merge_width)) {
      const std::size_t middle = std::min(left + merge_width, size);
      const std::size_t right = std::min(left + (2 * merge_width), size);
```

## 5. Детали реализации

`ValidationImpl` проверяет непустой вход. `PreProcessingImpl` копирует вход в
`data_` и очищает выход. `RunImpl` сортирует блоки и выполняет итеративное
слияние. `PostProcessingImpl` проверяет `data_` через
`std::ranges::is_sorted` и только после этого записывает `GetOutput()`.

При сортировке блок `block_start..block_end` не пересекается с соседним, потому
что шаг равен `kBlockSize`. При слиянии каждый поток пишет в свой диапазон
`merged_data[left, right)`. Общий `data` на фазе слияния только читается, а
`data.swap(merged_data)` выполняется после завершения параллельного цикла.

## 6. Проверка корректности

Эталон в функциональных тестах строится через `std::ranges::sort`. Запуск
`PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_func_tests` с фильтром
`*olesnitskiy_v_hoare_sort_simple_merge*` прошел: 75 tests, 60 passed, 15
ALL-tests skipped вне `mpirun`. ALL отдельно прошел под `mpirun -np 2`: 15
passed.

## 7. Экспериментальная среда

- **Сборка:** `build_olesnitskiy`
- **Compiler:** `g++-14`
- **Flags:** `-O3 -DNDEBUG`, `std=gnu++23`
- **Размер входных данных:** `N=100000`
- **Baseline TaskRun:** `0.0058254364 s`
- **Baseline pipeline:** `0.0068995056 s`
- **Число повторов:** 5 по умолчанию

Потоки задавались через `PPC_NUM_THREADS`.

## 8. Результаты

- threads: 1; time: 0.2042550788 s; speedup: 0.029; efficiency: 0.029; notes:
  `TaskRun`, `N=100000`.
- threads: 2; time: 0.2022266690 s; speedup: 0.029; efficiency: 0.014; notes:
  `TaskRun`, `N=100000`.
- threads: 4; time: 0.1880041056 s; speedup: 0.031; efficiency: 0.008; notes:
  `TaskRun`, `N=100000`.
- threads: 1; time: 0.1145523454 s; speedup: 0.060; efficiency: 0.060; notes:
  `pipeline`, `N=100000`.
- threads: 4; time: 0.1193738506 s; speedup: 0.058; efficiency: 0.014; notes:
  `pipeline`, `N=100000`.

## 9. Выводы

На данном размере входа OpenMP-версия медленнее baseline: лучший `TaskRun`
замер `0.1880041056 s` при 4 потоках дает speedup `0.031`. Основная причина —
много мелких блоков `kBlockSize=64` и дополнительные временные векторы при
слиянии.
