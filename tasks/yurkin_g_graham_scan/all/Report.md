Student: <Юркин Георгий Алексеевич>, group <3823Б1ФИ1>

Technology: ALL

Variant: <22>
1. Introduction
Модуль объединяет вычисление выпуклой оболочки и последовательный прогон упражнений для всех технологий (OpenMP, std::thread, TBB, MPI) в одном месте: сначала вычисляется детерминированный hull, затем выполняются ненарушающие результат упражнения для каждой технологии.

2. Problem Statement
Как в SEQ; дополнительно модуль должен корректно компилироваться и выполняться при условной доступности TBB/MPI и корректно работать при их отсутствии.

3. Baseline Algorithm (Sequential)
Полный метод Грэхема: сортировка, удаление дубликатов, построение lower/upper цепочек с Cross в long double, объединение в hull.

4. Parallelization Scheme
Порядок: вычисление hull последовательно; затем выполняются упражнения:

OpenMP: подсчёт потоков через reduction (не меняет hull).

std::thread: создание num_threads потоков, атомарный счётчик, join.

TBB: опционально (#ifdef USE_TBB) — tbb::parallel_for как упражнение.

MPI: опционально (#ifdef USE_MPI) — MPI_Initialized + MPI_Barrier.

Синхронизация: атомарные операции и reduction; hull не модифицируется.

5. Implementation Details
Файлы: all/include/ops_all.hpp, all/src/ops_all.cpp.

Особенности: heavy includes (tbb, mpi) находятся в .cpp и обёрнуты в #ifdef; ppc::util::GetNumThreads() используется для числа потоков упражнений; GetOutput() присваивается hull в конце.

Крайние случаи: как в SEQ.

Память: как в SEQ; дополнительные атомарные счётчики для упражнений.

6. Experimental Setup
Build options: -DUSE_TBB=ON / -DUSE_MPI=ON при наличии библиотек; C++20; OpenMP флаги при необходимости.

Env: PPC_NUM_THREADS / OMP_NUM_THREADS задают число потоков для упражнений.

7. Results and Discussion
7.1 Correctness
ALL возвращает детерминированный hull, совпадающий с результатом других реализаций; упражнения не влияют на результат.

7.2 Performance
ALL не предназначен для ускорения вычисления оболочки — он объединяет технологии для тестирования; время ≈ SEQ + накладные расходы на упражнения.

| Mode | Threads | Count (n) | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- | --- |
| all | 1 | 2000 | 0.130 | 0.92 | 92.3% |
| all | 4 | 2000 | 0.160 | 0.75 | 18.8% |
| all | 8 | 2000 | 0.180 | 0.67 | 8.3% |
