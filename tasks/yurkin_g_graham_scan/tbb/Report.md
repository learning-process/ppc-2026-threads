Student: <Юркин Георгий Алексеевич>, group <3823Б1ФИ1>

Technology: TBB

Variant: <22>

# TBB

## 1. Introduction

Цель — использовать Intel TBB для распараллеливания подготовительных этапов
(в первую очередь сортировки) и получить ускорение на больших входах.

## 2. Problem Statement

Как в SEQ.

## 3. Baseline Algorithm (Sequential)

Алгоритм тот же; TBB применяется для параллельной сортировки и подготовки,
финальный проход остаётся последовательным.

## 4. Parallelization Scheme

Декомпозиция: tbb::parallel_sort для сортировки в PreProcessingImpl и RunImpl;
tbb::parallel_for — опционально для вспомогательных вычислений.

Синхронизация: атомарные счётчики для упражнений; финальный проход —
последовательный.

## 5. Implementation Details

Файлы: tbb/include/ops_tbb.hpp, tbb/src/ops_tbb.cpp.

Особенности: tbb::parallel_sort используется для ускорения сортировки;
удаление дубликатов — последовательный проход; Cross в long double.

Память: TBB может использовать дополнительные временные буферы при сортировке.

## 6. Experimental Setup

- CPU: AMD Ryzen 5 5500U
- RAM: 16 ГБ
- OS: Windows 10

Toolchain: gcc/clang с C++20, Release (-O3).

Build: -DUSE_TBB=ON, find_package(TBB), линковка TBB::tbb.

Env: контроль числа потоков через TBB‑настройки или переменные среды.

## 7. Results and Discussion

### 7.1 Correctness

Результат совпадает с эталоном; параллельная сортировка не меняет семантику.

### 7.2 Performance

TBB даёт выигрыш на больших n в этапах сортировки; финальный проход остаётся
узким местом.

| Mode | Threads | Count (n) | Time, s | Speedup | Efficiency |
|------|---------|-----------|---------|---------|------------|
| tbb  | 1       | 2000      | 0.119   | 1.01    | 100.8%     |
| tbb  | 2       | 2000      | 0.082   | 1.46    | 73.2%      |
| tbb  | 4       | 2000      | 0.070   | 1.71    | 42.9%      |
| tbb  | 8       | 2000      | 0.058   | 2.07    | 25.9%      |

## 8. Conclusions

TBB эффективен при больших объёмах данных; использовать при наличии библиотеки
в CI.

## 9. References

Intel TBB документация.
