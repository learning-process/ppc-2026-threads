# Выделение рёбер на изображении с использованием оператора Собеля - ALL

- Student: Крюков Е.
- Technology: MPI, OMP
- Variant: 27

## 1. Контекст

Гибридная версия сочетает два уровня параллелизма: MPI распределяет работу
между процессами, OpenMP параллелит вычисления внутри каждого процесса.
Это позволяет эффективнее использовать ресурсы многоядерных систем.
Baseline - `seq/report.md`.

## 2. Постановка задачи

Идентична SEQ-версии: `Image{width, height, data}` - `std::vector<int>`.
Корректность проверяется совпадением с SEQ на всём тестовом наборе.

## 3. Базовый алгоритм

Подробно описан в `seq/report.md`. RGB - grayscale - свёртка ядром Собеля 3×3
по внутренним пикселям - `magnitude = sqrt(Gx² + Gy²)`.

## 4. Межпроцессная схема

Вся работа с MPI сосредоточена в `RunImpl`.
`ValidationImpl`, `PreProcessingImpl`, `PostProcessingImpl` идентичны SEQ -
каждый процесс самостоятельно готовит полный массив `grayscale_`
и выходной буфер.

**Роли rank-ов:** все процессы равноправны - нет выделенного мастер-процесса
для вычислений. Каждый rank вычисляет свой диапазон строк независимо,
без межпроцессных обменов в ходе вычислений.

**Распределение строк:** внутренние строки `[1, h-2]` делятся между `nproc`
процессами:

```cpp
// File: all/src/ops_all.cpp - RunImpl
int rank, nproc;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &nproc);

const int base  = (h - 2) / nproc;
const int extra = (h - 2) % nproc;

const int local_count = base + (rank < extra ? 1 : 0);
const int local_start = 1 + base * rank + std::min(rank, extra);
```

**MPI-вызовы и их назначение:**

| Вызов             | Место   | Назначение                                      |
| ----------------- | ------- | ----------------------------------------------- |
| `MPI_Comm_rank`   | RunImpl | Определить номер текущего процесса              |
| `MPI_Comm_size`   | RunImpl | Определить общее число процессов                |
| `MPI_Gatherv`     | RunImpl | Собрать локальные результаты на rank 0          |
| `MPI_Bcast`       | RunImpl | Разослать полный результат всем rank-ам         |

`MPI_Gatherv` использует `counts[i] = local_count_i * w` и
`displs[i] = local_start_i * w`, что позволяет записать локальные результаты
сразу в правильные позиции общего буфера `output`.

`MPI_Bcast` после `Gatherv` нужен, чтобы все процессы держали актуальный
`output` - это требование фреймворка.

## 5. Внутрипроцессная схема

Внутри каждого rank-а цикл по локальным строкам параллелится через OpenMP:

```cpp
// File: all/src/ops_all.cpp - RunImpl (OpenMP-фрагмент)
#pragma omp parallel for default(none) \
    shared(local_output, gray, gx_kernel, gy_kernel) \
    firstprivate(local_start, local_count, w) \
    schedule(static)
for (int li = 0; li < local_count; ++li) {
  const int row = local_start + li;
  for (int col = 1; col < w - 1; ++col) {
    // вычисление gx, gy, magnitude
    local_output[li * w + col] = magnitude;
  }
}
```

Атрибуты переменных аналогичны OMP-версии. Неявный барьер OpenMP в конце
`parallel for` гарантирует заполнение `local_output` до вызова `MPI_Gatherv`.

**Конфигурация:** 4 MPI-процесса x потоки по `OMP_NUM_THREADS` (4 в тестах)
= 16 рабочих единиц суммарно.

## 6. Детали реализации

Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`.

Заголовочный файл идентичен SEQ/OMP по составу полей
(`grayscale_`, `width_`, `height_`). MPI-специфичные переменные
(`rank`, `nproc`, `counts`, `displs`, `local_output`) объявлены локально
в `RunImpl` - это изолирует логику MPI от остальных методов.

**Потенциальные узкие места:**

- `MPI_Gatherv`: пересылка `w * local_count * sizeof(int)` байт от каждого
  rank-а на rank 0. При `size=2048`, `nproc=4` - около 4 МБ суммарно.
- `MPI_Bcast`: рассылка полного выходного буфера (`w*h*sizeof(int)`) всем
  процессам. При `size=2048` - около 16 МБ.

## 7. Проверка корректности

Все три функциональных теста (ConstantImage, VerticalEdge, HorizontalEdge)
пройдены при запуске с `mpiexec -n 4`.
Результат побайтово совпадает с SEQ.
При `nproc=1` поведение сводится к однопроцессному OMP - результат не меняется.

Согласованность между rank-ами обеспечивается `MPI_Bcast`: после его завершения
`GetOutput()` идентичен на всех процессах.

## 8. Экспериментальная среда

| Параметр                       | Значение                             |
| ------------------------------ | ------------------------------------ |
| CPU                            | AMD Ryzen 3 3100 (4 физических ядра) |
| ОС                             | Windows 10                           |
| Компилятор                     | MSVC (Release)                       |
| CMake build type               | Release                              |
| `PPC_NUM_PROC` / `mpiexec -n`  | 4                                    |
| `OMP_NUM_THREADS`              | 4                                    |
| Конфигурация                   | 4 ranks x 4 OMP-потока = 16 единиц  |

Команда запуска:

```bash
mpiexec -n 4 .\build\bin\ppc_perf_tests.exe --gtest_filter="*KrykovE*"
```

## 9. Результаты

Ускорение считается по `T_seq_task_run / T_all_task_run`.
Эффективность `E = speedup / total_workers`,
где `total_workers = ranks x threads_per_rank = 16`.

| size | SEQ task_run, с | ALL task_run, с | Speedup | Workers | Efficiency |
| ---- | --------------- | --------------- | ------- | ------- | ---------- |
| 512  | 0.00738         | 0.00230         | 3.21x   | 16      | 20.1%      |
| 1024 | 0.03471         | 0.01281         | 2.71x   | 16      | 16.9%      |
| 2048 | 0.11435         | 0.03637         | 3.14x   | 16      | 19.6%      |

| size | SEQ pipeline, с | ALL pipeline, с | Speedup |
| ---- | --------------- | --------------- | ------- |
| 512  | 0.00851         | 0.00379         | 2.24x   |
| 1024 | 0.03374         | 0.01320         | 2.56x   |
| 2048 | 0.12496         | 0.05271         | 2.37x   |

ALL-версия показывает наилучшее абсолютное ускорение среди всех технологий:
до 3.21x по `task_run`. Разрыв между `pipeline` и `task_run` объясняется тем,
что `pipeline`-режим включает последовательный RGB - grayscale и коммуникации
MPI в overhead.

## 10. Выводы

Гибридная схема MPI + OpenMP оправдана при больших изображениях (>=512x512):
накладные расходы на MPI-коммуникации компенсируются реальным параллелизмом
двух уровней. На данной задаче ALL показывает лучшее абсолютное время из всех
пяти технологий. Основные ограничения - memory-bound характер алгоритма
и линейный рост коммуникационного overhead с размером изображения.
