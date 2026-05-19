# Отчет OMP: сортировка Хоара с простым слиянием

## Контекст и базовый алгоритм

Вход и выход: `std::vector<int>`
(`common/include/common.hpp:11`).
Последовательное ядро: разбиение Хоара и quicksort по локальным диапазонам
(`omp/src/ops_omp.cpp:21`). После локальной сортировки
блоков по 64 элемента выполняется попарное простое слияние до полного массива
(`omp/src/ops_omp.cpp:115`). Средняя сложность локальной
quicksort-части `O(n log n)`, слияния по уровням дают дополнительный линейный
проход на каждом уровне.

## Схема распараллеливания

OpenMP используется в двух `parallel for`: сортировка независимых блоков и
слияние независимых пар блоков (`omp/src/ops_omp.cpp:119`,
`omp/src/ops_omp.cpp:130`). В директивах стоит
`default(none)`, поэтому все shared-переменные перечислены явно: `size`, `data`,
`merged_data`, `merge_width`. Индексы циклов и временные векторы являются
private по области видимости тела цикла. Reduction не нужен: нет общей скалярной
агрегации. Явный `schedule` не задан, значит применяется runtime-default OpenMP;
для равномерных блоков по 64 элемента это соответствует идее статического
распределения без дополнительных вычислений планировщика. Барьер OpenMP в конце
каждого `parallel for` нужен перед `data.swap(merged_data)`.

Фрагмент, `omp/src/ops_omp.cpp:119`: параллельные блоки
пишут в непересекающиеся диапазоны.

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

## Детали pipeline

`ValidationImpl` проверяет непустой вход
(`omp/src/ops_omp.cpp:100`). `PreProcessingImpl` копирует
вход в `data_` и очищает выход (`omp/src/ops_omp.cpp:104`).
`RunImpl` сортирует блоки и выполняет итеративное слияние
(`omp/src/ops_omp.cpp:110`). `PostProcessingImpl` проверяет
`data_` через `std::ranges::is_sorted` и только после этого записывает
`GetOutput()` (`omp/src/ops_omp.cpp:156`).

## Проверка отсутствия гонок

При сортировке блок `block_start..block_end` не пересекается с соседним, потому
что шаг равен `kBlockSize` (`omp/src/ops_omp.cpp:120`). При
слиянии каждый поток пишет в `merged_data[left, right)`, где `left` идет с шагом
`2 * merge_width` (`omp/src/ops_omp.cpp:131`). Общий `data`
на фазе слияния только читается, а `data.swap(merged_data)` выполняется после
завершения параллельного цикла (`omp/src/ops_omp.cpp:150`).

## Корректность и среда

Эталон в функциональных тестах строится `std::ranges::sort`
(`tests/functional/main.cpp:35`). Запуск
`PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_func_tests` с фильтром
`*olesnitskiy_v_hoare_sort_simple_merge*` прошел: 75 tests, 60 passed, 15
ALL-tests skipped вне `mpirun`. ALL отдельно прошел под `mpirun -np 2`: 15
passed. Performance-вход: 100000 случайных чисел
(`tests/performance/main.cpp:20`).

## Результаты

Baseline: `seq` `TaskRun = 0.0058254364 s`; для pipeline baseline `0.0068995056
s`. Потоки задавались через `PPC_NUM_THREADS`; раннер ограничивает TBB этой же
переменной, а OpenMP использует ее через окружение сборки/запуска. Framework
выполняет 5 повторов по умолчанию
(`modules/performance/include/performance.hpp:21`).

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

## Выводы

На данном размере входа OpenMP-версия медленнее baseline: лучший замер
`0.1880041056 s` при 4 потоках дает speedup `0.031`. Причина подтверждается
кодом: блок `kBlockSize=64` создает много мелких задач и временных векторов при
слиянии (`omp/src/ops_omp.cpp:115`,
`omp/src/ops_omp.cpp:136`).
