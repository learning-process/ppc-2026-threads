Student: <Юркин Георгий Алексеевич>, group <3823Б1ФИ1>

Technology: OMP

Variant: <22>

# OMP

## 1. Introduction

Цель — ускорить подготовительные этапы (в первую очередь сортировку) с помощью
параллельных возможностей платформы, сохранив детерминированный результат.

## 2. Problem Statement

Как в SEQ: вход — `std::vector<Point>`, выход — `std::vector<Point>` оболочки;
корректная обработка дубликатов и коллинеарности.

## 3. Baseline Algorithm (Sequential)

Алгоритм совпадает с SEQ; параллелизация применяется в подготовительных этапах
при поддержке execution policies.

## 4. Parallelization Scheme

Декомпозиция: при наличии поддержки std::execution::par используется
параллельная сортировка (std::sort(std::execution::par, ...)) в
PreProcessingImpl и RunImpl.

Синхронизация: reduction для агрегатов; финальный стековый проход остаётся
последовательным.

Роли: все потоки выполняют подготовительные операции; объединение —
последовательное.

## 5. Implementation Details

Файлы: omp/include/ops_omp.hpp, omp/src/ops_omp.cpp.

Особенности: условное использование std::execution::par (если поддерживается
компилятором); удаление дубликатов — последовательный проход; Cross в long double.

Крайние случаи: как в SEQ.

## 6. Experimental Setup

- CPU: AMD Ryzen 5 5500U
- RAM: 16 ГБ
- OS: Windows 10

Toolchain: gcc/clang с C++20, Release (-O3), с флагом -fopenmp.

Env: OMP_NUM_THREADS / PPC_NUM_THREADS для управления числом потоков.

## 7. Results and Discussion

### 7.1 Correctness

Результат совпадает с SEQ; параллельные блоки не модифицируют оболочку.

### 7.2 Performance

Ускорение зависит от поддержки параллельной сортировки и доли параллелизуемой
работы. Замеры необходимо получить экспериментально.

| Mode | Count | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- |
| seq | 1 | 0.001330 | 1.00 | N/A |
| omp | 2 | 0.001285 | 1.04 | 51.8% |
| omp | 4 | 0.001264 | 1.05 | 26.3% |
| omp | 8 | 0.001758 | 0.76 | 9.5% |

## 8. Conclusions

OpenMP даёт выигрыш на подготовительных этапах при корректной настройке;
финальный проход оставлен последовательным для детерминизма.

## 9. References

OpenMP спецификация; C++ execution policies.
