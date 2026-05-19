# Отчет SEQ: сортировка Хоара с простым слиянием

## Контекст и постановка

Задача: отсортировать `std::vector<int>` по неубыванию. Типы входа и выхода
заданы как `std::vector<int>` в
`common/include/common.hpp:11`.
Последовательная версия использует разбиение Хоара и итеративный quicksort:
опорный элемент берется из середины диапазона, индексы двигаются навстречу,
элементы меняются местами до пересечения
(`seq/src/ops_seq.cpp:18`). Средняя сложность сортировки:
`O(n log n)`, худшая: `O(n^2)`, память стека диапазонов: до `O(n)` в
неблагоприятном случае. Инвариант разбиения: после возврата `j` элементы слева
от границы не больше опорной группы, справа не меньше; рекурсивные подзадачи
заменены стеком диапазонов (`seq/src/ops_seq.cpp:42`).

Фрагмент реализации, `seq/src/ops_seq.cpp:18`: выбор pivot и
схема Хоара.

```cpp
int OlesnitskiyVHoareSortSimpleMergeSEQ::HoarePartition(
    std::vector<int> &values, int left, int right) {
  const int pivot = values[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    ++i;
    while (values[i] < pivot) {
      ++i;
    }

    --j;
    while (values[j] > pivot) {
      --j;
    }

    if (i >= j) {
      return j;
    }

    std::swap(values[i], values[j]);
  }
}
```

## Детали реализации

`ValidationImpl` принимает только непустой вход
(`seq/src/ops_seq.cpp:66`). `PreProcessingImpl` копирует
вход в выходной буфер, чтобы сортировать результат без изменения `GetInput()`
(`seq/src/ops_seq.cpp:70`). `RunImpl` пропускает массивы
размера `0/1`, затем вызывает `HoareQuickSort` и проверяет
`std::ranges::is_sorted` (`seq/src/ops_seq.cpp:75`).
`PostProcessingImpl` повторно проверяет непустой и отсортированный выход
(`seq/src/ops_seq.cpp:86`).

## Проверка корректности

Функциональный тест строит эталон через `std::ranges::sort(expected_data)` и
сравнивает его с выходом задачи
(`tests/functional/main.cpp:35`). Набор
содержит 15 сценариев: один элемент, дубликаты, обратный порядок, отрицательные
значения, границу блока 64 и `int`-пределы
(`tests/functional/main.cpp:59`). Запуск
`PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_func_tests` с фильтром
`*olesnitskiy_v_hoare_sort_simple_merge*` прошел: 75 тестов, 60 passed и 15
ALL-тестов skipped вне `mpirun`. Отдельный запуск `mpirun -np 2
./build_olesnitskiy/bin/ppc_func_tests` с фильтром
`*olesnitskiy_v_hoare_sort_simple_merge_all_enabled*` прошел 15/15 ALL-тестов.

## Экспериментальная среда

Сборка `build_olesnitskiy`: `g++-14`, `-O3 -DNDEBUG`, `std=gnu++23`, OpenMPI и
oneTBB видны в
`build_olesnitskiy/compile_commands.json`.
Размер performance-набора: 100000 случайных `int` из `[-1000000, 1000000]`
(`tests/performance/main.cpp:20`). Время
печатает `Perf::PrintPerfStatistic` после `TaskRun` или `PipelineRun`
(`modules/util/include/perf_test_util.hpp:87`).
Число повторов по умолчанию равно 5
(`modules/performance/include/performance.hpp:21`).
Важно для интерпретации: `TaskRun` повторяет только `Run()` после одного
`PreProcessing()`
(`modules/performance/include/performance.hpp:62`),
а вход для каждого gtest-параметра генерируется заново через
`std::random_device`
(`tests/performance/main.cpp:24`).

## Результаты

Baseline измерен командой: `PPC_NUM_THREADS=1
./build_olesnitskiy/bin/ppc_perf_tests` с фильтром
`*task_run*olesnitskiy_v_hoare_sort_simple_merge_seq_enabled`.

- backend: seq; time: 0.0058254364 s; speedup: 1.000; efficiency: 1.000; notes:
  `N=100000`, `TaskRun`, `PPC_NUM_THREADS=1`.
- backend: seq; time: 0.0068995056 s; speedup: 1.000; efficiency: 1.000; notes:
  `N=100000`, `pipeline`, `PPC_NUM_THREADS=1`.

## Выводы

Последовательная версия является baseline для speedup. Корректность подтверждена
сравнением с `std::ranges::sort`; производительность подтверждена `TaskRun` и
`pipeline`-замерами на `N=100000`.
