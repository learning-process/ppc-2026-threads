Student: <Юркин Георгий Алексеевич>, group <3823Б1ФИ1>

Technology: STL

Variant: <22>

# STL

## 1. Introduction

Реализация использует стандартные средства C++20 (std::ranges) для компактного,
переносимого и детерминированного кода.

## 2. Problem Statement

Как в SEQ.

## 3. Baseline Algorithm (Sequential)

Полная реализация метода Грэхема: сортировка (std::ranges::sort),
удаление дубликатов, построение lower/upper цепочек с Cross в long double,
объединение в hull.

## 4. Parallelization Scheme

Не применяется — чисто последовательная STL‑реализация.

## 5. Implementation Details

Файлы: stl/include/ops_stl.hpp, stl/src/ops_stl.cpp.

Особенности: std::ranges::sort, std::ranges::reverse_view;
при несовместимости toolchain заменить на std::sort.

Крайние случаи: обработаны как в SEQ.

## 6. Experimental Setup

- CPU: AMD Ryzen 5 5500U
- RAM: 16 ГБ
- OS: Windows 10

Toolchain: gcc/clang с C++20, Release (-O3).

## 7. Results and Discussion

### 7.1 Correctness

Проходит функциональные тесты; корректно обрабатывает дубликаты и коллинеарность.

### 7.2 Performance

Время аналогично SEQ (эталон). Замеры — по результатам локального профилирования.

| Mode | Count | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- |
| seq | 1 | 0.001330 | 1.00 | N/A |
| stl | 2 | 0.001346 | 0.99 | 49.5% |
| stl | 4 | 0.001346 | 0.99 | 24.8% |
| stl | 8 | 0.002058 | 0.65 | 8.1% |

## 8. Conclusions

Чистая, читабельная и переносимая реализация; рекомендуется как эталон.

## 9. References

C++20 стандарт, документация по std::ranges.
