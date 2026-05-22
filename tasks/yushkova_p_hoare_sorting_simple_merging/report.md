# Сортировка Хоара с простым слиянием - сводный отчет

- **Студент:** Юшкова Полина Александровна, 3823Б1ПР2
- **Технологии:** SEQ, OMP, TBB, STL, ALL (MPI + STL)
- **Вариант:** 13

## 1. Контекст

Работа сравнивает реализации сортировки `std::vector<int>` по неубыванию для backend-ов `seq`, `omp`, `tbb`, `stl` и
`all`. Последовательная версия задает baseline, потоковые версии проверяют разные модели внутрипроцессного
распараллеливания, а ALL добавляет MPI-уровень.

## 2. Постановка задачи

- **Входные данные:** непустой объект `std::vector<int>`.
- **Выходные данные:** объект `std::vector<int>` того же размера, элементы которого переставлены в порядке неубывания.
- **Эталон корректности:** совпадение с результатом `std::ranges::sort` в функциональных тестах.
- **Baseline:** `seq`, `T_seq(task_run) = 0.0030497600 s`, `T_seq(pipeline) = 0.0083235000 s`.

Единые типы входа и выхода заданы в `common/include/common.hpp`.

## 3. Базовый алгоритм

Базовое вычислительное ядро использует сортировку Хоара: опорный элемент берется из середины диапазона, затем диапазон
разбивается двумя индексами и сортируется итеративным quicksort.

Параллельные реализации используют принцип простого слияния:

- локальная сортировка независимых частей (половин или блоков по 64 элемента);
- уровни слияния объединяют соседние отсортированные диапазоны.

## 4. Схемы распараллеливания

- **SEQ:** две половины, сортировка Хоара, `SimpleMerge`.
- **OMP:** деление массива на `chunks`, сортировка блоков в `#pragma omp parallel for`, затем последовательное слияние.
- **TBB:** `oneapi::tbb::parallel_for(blocked_range<size_t>)` для сортировки блоков по 64 и для каждого прохода слияния.
- **STL:** при `PPC_NUM_THREADS > 1` используется 2 worker-а (один `std::thread` + текущий поток) для сортировки двух
  половин, затем merge.
- **ALL:** MPI распределяет фрагменты между rank-ами (`Scatterv/Gatherv/Bcast`), а локальная часть каждого rank-а
  сортирует фрагмент STL-схемой (блоки по 64 и уровни слияния).

## 5. Детали методики

Замеры выполнены на `N=100000` случайных `int` из диапазона `[-1000000, 1000000]`, число повторов - 5.

- `TaskRun`: повторяется только `Run()` после одного `PreProcessing()`.
- `pipeline`: каждый повтор проходит полный pipeline (`Validation` → `PreProcessing` → `Run` → `PostProcessing`).

Speedup считался как `T_seq / T_backend`, efficiency - как `speedup / workers`.

## 6. Проверка корректности

Функциональные тесты регистрируют все backend-ы (`seq/omp/stl/tbb/all`) и сравнивают результат с `std::ranges::sort`.
Дополнительно в реализациях используются проверки `std::ranges::is_sorted`.

## 7. Экспериментальная среда

- **ОС:** Windows
- **Размер входных данных:** `N=100000`
- **Диапазон значений:** `[-1000000, 1000000]`
- **Число повторов:** 5
- **ALL:** `mpiexec` (Microsoft MPI)
- **Auto workers (ALL локально):** `std::thread::hardware_concurrency() = 12` на тестовой машине

## 8. Результаты TaskRun

- backend: seq; time: 0.0030497600 s; speedup: 1.000; efficiency: 1.000; notes: baseline.
- backend: omp, 1 thread; time: 0.0029892200 s; speedup: 1.020; efficiency: 1.020; notes: `PPC_NUM_THREADS=1`.
- backend: omp, 2 threads; time: 0.0021860800 s; speedup: 1.395; efficiency: 0.698; notes: `PPC_NUM_THREADS=2`.
- backend: omp, 4 threads; time: 0.0018365600 s; speedup: 1.661; efficiency: 0.415; notes: `PPC_NUM_THREADS=4`.
- backend: stl, workers: 1; time: 0.0028924000 s; speedup: 1.054; efficiency: 1.054; notes: `PPC_NUM_THREADS=1`.
- backend: stl, workers: 2; time: 0.0021834000 s; speedup: 1.397; efficiency: 0.698;
  notes: `PPC_NUM_THREADS=4` включает параллельную ветку.
- backend: tbb, 1 worker; time: 0.0076817000 s; speedup: 0.397; efficiency: 0.397;
  notes: `PPC_NUM_THREADS=1`, `global_control`.
- backend: tbb, 2 workers; time: 0.0048485400 s; speedup: 0.629; efficiency: 0.315;
  notes: `PPC_NUM_THREADS=2`, `global_control`.
- backend: tbb, 4 workers; time: 0.0028551400 s; speedup: 1.068; efficiency: 0.267;
  notes: `PPC_NUM_THREADS=4`, `global_control`.
- backend: all; ranks: 1; threads_per_rank: 12; total_workers: 12; time: 0.0159865600 s;
  speedup: 0.191; efficiency: 0.016; notes: `mpiexec -n 1`.
- backend: all; ranks: 2; threads_per_rank: 12; total_workers: 24; time: 0.0155627800 s;
  speedup: 0.196; efficiency: 0.008; notes: `mpiexec -n 2`.
- backend: all; ranks: 4; threads_per_rank: 12; total_workers: 48; time: 0.0125462200 s;
  speedup: 0.243; efficiency: 0.005; notes: `mpiexec -n 4`.

## 9. Результаты Pipeline

- backend: seq; time: 0.0083235000 s; speedup: 1.000; efficiency: 1.000; notes: baseline.
- backend: omp, 1 thread; time: 0.0112251400 s; speedup: 0.742; efficiency: 0.742; notes: `PPC_NUM_THREADS=1`.
- backend: omp, 4 threads; time: 0.0037493600 s; speedup: 2.220; efficiency: 0.555; notes: `PPC_NUM_THREADS=4`.
- backend: stl, workers: 1; time: 0.0086819200 s; speedup: 0.959; efficiency: 0.959; notes: `PPC_NUM_THREADS=1`.
- backend: stl, workers: 2; time: 0.0055119600 s; speedup: 1.510; efficiency: 0.755;
  notes: `PPC_NUM_THREADS=4`.
- backend: tbb, 1 worker; time: 0.0113887600 s; speedup: 0.731; efficiency: 0.731;
  notes: `PPC_NUM_THREADS=1`, `global_control`.
- backend: tbb, 4 workers; time: 0.0037098400 s; speedup: 2.244; efficiency: 0.561;
  notes: `PPC_NUM_THREADS=4`, `global_control`.
- backend: all; ranks: 1; threads_per_rank: 12; total_workers: 12; time: 0.0173728400 s;
  speedup: 0.479; efficiency: 0.040; notes: `mpiexec -n 1`.
- backend: all; ranks: 4; threads_per_rank: 12; total_workers: 48; time: 0.0142562000 s;
  speedup: 0.584; efficiency: 0.012; notes: `mpiexec -n 4`.

## 10. Интерпретация результатов

`seq` служит baseline. `omp` показывает ускорение на `TaskRun` и заметный выигрыш на pipeline при 4 потоках, так как часть
работы распараллеливается по независимым блокам. `stl` ускоряет сортировку половин только при включении параллельной
ветки (2 worker-а). `tbb` на `N=100000` чувствителен к накладным расходам `parallel_for` и выделения буферов слияния:
в `TaskRun` speedup заметен только при 4 потоках, а в pipeline - уверенный выигрыш при 4 потоках.
`all` добавляет стоимость `Scatterv/Gatherv/Bcast` и финальное слияние на rank 0, из-за чего эффективность по
`total_workers` остается низкой.

## 11. Выводы

Для измеренного размера `N=100000` лучший pipeline-результат среди потоковых backend-ов показали `tbb(4)` и `omp(4)`
(времена порядка `0.0037 s`). ALL-версия на `4` rank-ах не ускоряет baseline на этом размере из-за MPI-обменов и
накладных расходов локальной сортировки.
