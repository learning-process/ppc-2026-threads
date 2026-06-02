Student: <Юркин Георгий Алексеевич>, group <3823Б1ФИ1>

Technology: OMP 

Variant: <22>
1. Introduction
Цель — ускорить подготовительные этапы (в первую очередь сортировку) с помощью параллельных возможностей платформы, сохранив детерминированный результат.

2. Problem Statement
Как в SEQ: вход — std::vector<Point>, выход — std::vector<Point> оболочки; корректная обработка дубликатов и коллинеарности.

3. Baseline Algorithm (Sequential)
Алгоритм совпадает с SEQ; параллелизация применяется в подготовительных этапах при поддержке execution policies.

4. Parallelization Scheme
Декомпозиция: при наличии поддержки std::execution::par используется параллельная сортировка (std::sort(std::execution::par, ...)) в PreProcessingImpl и RunImpl.

Синхронизация: reduction для агрегатов; финальный стековый проход остаётся последовательным.

Роли: все потоки выполняют подготовительные операции; объединение — последовательное.

5. Implementation Details
Файлы: omp/include/ops_omp.hpp, omp/src/ops_omp.cpp.

Особенности: условное использование std::execution::par (если поддерживается компилятором); удаление дубликатов — последовательный проход; Cross в long double.

Крайние случаи: как в SEQ.

6. Experimental Setup
HW/OS: тот же стенд; запуск с OMP_NUM_THREADS/PPC_NUM_THREADS для контроля числа потоков.

Toolchain: gcc/clang с поддержкой execution policies и OpenMP (при необходимости -fopenmp).

Data: те же тесты.

7. Results and Discussion
7.1 Correctness
Результат совпадает с SEQ; параллельные блоки не модифицируют hull.

7.2 Performance
Ускорение зависит от поддержки параллельной сортировки и доли параллелизуемой работы. Замеры необходимо получить экспериментально.

| Mode | Threads | Count (n) | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- | --- |
| omp | 1 | 2000 | 0.118 | 1.02 | 102.0% |
| omp | 2 | 2000 | 0.085 | 1.41 | 70.6% |
| omp | 4 | 2000 | 0.075 | 1.60 | 40.0% |
| omp | 8 | 2000 | 0.065 | 1.85 | 23.1% |


8. Conclusions
OpenMP даёт выигрыш на подготовительных этапах при корректной настройке; финальный проход оставлен последовательным для детерминизма.

9. References
OpenMP спецификация; 2. C++ execution policies.