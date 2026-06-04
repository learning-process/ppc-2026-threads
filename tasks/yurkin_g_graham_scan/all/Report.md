Student: <Юркин Георгий Алексеевич>, group <3823Б1ФИ1>

Technology: ALL

Variant: <22>

# ALL

## 1. Introduction

Модуль объединяет вычисление выпуклой оболочки и последовательный прогон
упражнений для всех технологий (OpenMP, std::thread, TBB, MPI) в одном месте:
сначала вычисляется детерминированный hull, затем выполняются ненарушающие
результат упражнения для каждой технологии.

## 2. Problem Statement

Как в SEQ; дополнительно модуль должен корректно компилироваться и выполняться
при условной доступности TBB/MPI и корректно работать при их отсутствии.

## 3. Baseline Algorithm (Sequential)

Полный метод Грэхема: сортировка, удаление дубликатов, построение lower/upper
цепочек с Cross в long double, объединение в hull.

## 4. Parallelization Scheme

Порядок: вычисление hull последовательно; затем выполняются упражнения:

OpenMP: подсчёт потоков через reduction (не меняет hull).

std::thread: создание num_threads потоков, атомарный счётчик, join.

TBB: опционально (#ifdef USE_TBB) — tbb::parallel_for как упражнение.

MPI: опционально (#ifdef USE_MPI) — MPI_Initialized + MPI_Barrier.

Синхронизация: атомарные операции и reduction; hull не модифицируется.

## 5. Implementation Details

Файлы: all/include/ops_all.hpp, all/src/ops_all.cpp.

Особенности: heavy includes (tbb, mpi) находятся в .cpp и обёрнуты в #ifdef;
ppc::util::GetNumThreads() используется для числа потоков упражнений;
GetOutput() присваивается hull в конце.

Крайние случаи: как в SEQ.

Память: как в SEQ; дополнительные атомарные счётчики для упражнений.

## 6. Experimental Setup

- CPU: AMD Ryzen 5 5500U
- RAM: 16 ГБ
- OS: Windows 10

Build options: -DUSE_TBB=ON / -DUSE_MPI=ON при наличии библиотек; C++20;

## 7. Results and Discussion

### 7.1 Correctness

ALL возвращает детерминированный hull, совпадающий с результатом других
реализаций; упражнения не влияют на результат.

### 7.2 Performance

ALL не предназначен для ускорения вычисления оболочки — он объединяет технологии
для тестирования; время ≈ SEQ + накладные расходы на упражнения.

| Mode | Count | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- |
| seq | 1 | 0.001330 | 1.00 | N/A |
| all | 2 | 0.004052 | 0.33 | 16.5% |
| all | 4 | 0.003923 | 0.34 | 8.5% |
| all | 8 | 0.005525 | 0.24 | 3.0% |

## 8. Conclusions

ALL полезен для интеграционного тестирования и демонстрации поддержки
технологий, но не оптимизирует сам алгоритм построения оболочки.

## 9. References

Intel TBB documentation.
OpenMP specification.
MPI standard.
