Student: <Юркин Георгий Алексеевич>, group <3823Б1ФИ1>

Technology: SEQ

Variant: <22>

# SEQ

## 1. Introduction

Построение выпуклой оболочки множества точек в 2D. Цель — получить
детерминированный вектор вершин оболочки в порядке обхода, корректно
обрабатывать дубликаты и коллинеарность.

## 2. Problem Statement

Вход: `std::vector<Point>` (Point { double x, double y }).
Выход: `std::vector<Point>` — вершины выпуклой оболочки в обходе.
Ограничения: произвольное количество точек; возможны дубликаты и коллинеарные
точки; корректная обработка случаев 0–2 точки.

## 3. Baseline Algorithm (Sequential)

Сортировка по x, затем y (std::ranges::sort).

Последовательное удаление точных дубликатов.

Если размер ≤ 2 — вернуть вход.

Построение lower и upper цепочек методом монотонной цепочки (Graham):
ориентация через Cross (в long double), условие удаления вершины
Cross(...) <= 0.

Объединение цепочек без дублирования крайних точек → hull.
Сложность: O(n log n). Память: O(n).

## 4. Parallelization Scheme

Не применяется (последовательная реализация).

## 5. Implementation Details

Файлы: seq/include/ops_seq.hpp, seq/src/ops_seq.cpp.

Ключевые функции: ValidationImpl, PreProcessingImpl, RunImpl,
PostProcessingImpl.

Вспомогательное: Cross в long double.

Крайние случаи: 0–2 точки, дубликаты, коллинеарность — обработаны.

Память: reserve для векторов lower, upper, hull.

## 6. Experimental Setup

- CPU: AMD Ryzen 5 5500U
- RAM: 16 ГБ
- OS: Windows 10

Toolchain: gcc/clang с C++20, Release (-O3).

Env: нет специальных переменных.

## 7. Results and Discussion

### 7.1 Correctness

Проверено unit‑тестами: фиксированные наборы (границы + внутренние точки)
возвращают ожидаемые вершины; PostProcessingImpl проверяет непустоту
результата при непустом входе.

### 7.2 Performance

| Mode | Count | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- |
| seq | 1 | 0.001330 | 1.00 | N/A |

## 8. Conclusions

Надёжная последовательная реализация, служащая базой для сравнения и проверки
корректности параллельных версий.

## 9. References

O'Rourke J., Computational Geometry in C.
C++20 standard (std::ranges).
