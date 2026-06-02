Student: <Юркин Георгий Алексеевич>, group <3823Б1ФИ1>

Technology: STL

Variant: <22>
1. Introduction
Реализация использует стандартные средства C++20 (std::ranges) для компактного, переносимого и детерминированного кода.

2. Problem Statement
Как в SEQ.

3. Baseline Algorithm (Sequential)
Полная реализация метода Грэхема: сортировка (std::ranges::sort), удаление дубликатов, построение lower/upper цепочек с Cross в long double, объединение в hull.

4. Parallelization Scheme
Не применяется — чисто последовательная STL‑реализация.

5. Implementation Details
Файлы: stl/include/ops_stl.hpp, stl/src/ops_stl.cpp.

Особенности: std::ranges::sort, std::ranges::reverse_view; при несовместимости toolchain заменить на std::sort.

Крайние случаи: обработаны как в SEQ.

6. Experimental Setup
Toolchain: требуется поддержка C++20 и std::ranges.

Data: те же тесты.

7. Results and Discussion
7.1 Correctness
Проходит функциональные тесты; корректно обрабатывает дубликаты и коллинеарность.

7.2 Performance
Время аналогично SEQ (эталон). Замеры — по результатам локального профилирования.

| Mode | Count (n) | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- |
| stl | 2000 | 0.118 | 1.02 | 101.7% |

8. Conclusions
Чистая, читабельная и переносимая реализация; рекомендуется как эталон.

9. References
C++20 стандарт, документация по std::ranges.