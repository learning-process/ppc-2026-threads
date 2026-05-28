# OMP: Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – столбцовый (CCS)

Обзор задачи, метрики и сборка: [../report.md](../report.md).

- Студент: Борунов Владислав Алексеевич
- Группа: 3823Б1ПР3
- Вариант: 7
- Задача: Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – столбцовый (CCS).
- Технологии: OMP

## 1. Назначение

BorunovVComplexCcsOMP - параллельное умножение через OpenMP: статическое разбиение столбцов B,
приватные аккумуляторы, слияние CCS после барьера parallel.

## 2. Распараллеливание

Столбцы C[:,j] независимы. Поток tid обрабатывает

```cpp
jstart = (tid * bc) / nt
jend   = ((tid + 1) * bc) / nt
```

bc = B.num_cols, nt = omp_get_num_threads().

Запись только в t_values[tid], t_row_indices[tid], t_col_nnz[tid] - конфликтов нет.

Столбец (ProcessColumn): обход B[:,j] и A[:,p] - acc[i] += ...; маркер marker[i] != j;
CustomShellSort(touched); порог |acc[i]| > 1e-9.

## 3. Синхронизация и OpenMP

```cpp
#pragma omp parallel default(none) \
    shared(a, b, bc, t_values, t_row_indices, t_col_nnz) \
    num_threads(ppc::util::GetNumThreads())
```

| Клауза           | Назначение                               |
|------------------|------------------------------------------|
| parallel         | команда потоков                          |
| default(none)    | запрет неявного sharing                  |
| shared(...)      | входные матрицы и векторы буферов по tid |
| num_threads(...) | GetNumThreads() - PPC_NUM_THREADS        |

Private (локально в теле): tid, nt, jstart, jend, acc, marker, touched.

Reduction: не используется; сумма в приватном acc, сброс в буфер потока.

Schedule: не задан; цикл по j распределён вручную (статическое разбиение по столбцам, без omp for).

Merge (MergeResults) - в главном потоке после parallel: сборка C.col_ptrs, reserve, insert из буферов.

## 4. Отличия от других реализаций

| Аспект            | OMP                  | SEQ            | TBB                          | STL                    |
|-------------------|----------------------|----------------|------------------------------|------------------------|
| workers           | PPC_NUM_THREADS      | 1              | PPC_NUM_THREADS              | hardware_concurrency() |
| Маркер            | marker[i]!=j         | `vector<bool>` | как OMP                      | `vector<bool>`         |
| Параллельный цикл | #pragma omp parallel | -              | parallel_for по tid          | std::thread            |
| Merge             | последовательный     | -              | col_ptrs + параллельный copy | после join()           |

## 5. Конфигурация workers

workers = ppc::util::GetNumThreads() (переменная окружения PPC_NUM_THREADS, иначе 1).

## 6. Результаты (локальный прогон)

Baseline: SEQ task_run = 0.0408020800 с., SEQ pipeline = 0.0332925600 с.

### 6.1 task_run

| Workers | time, с      | speedup | efficiency |
|---------|--------------|---------|------------|
| 4       | 0.0618384400 | 0.66    | 16.5%      |
| 8       | 0.0530478200 | 0.77    | 9.6%       |
| 16      | 0.0436728800 | 0.93    | 5.8%       |

### 6.2 pipeline

| Workers | time, с      | speedup | efficiency |
|---------|--------------|---------|------------|
| 4       | 0.0484773800 | 0.69    | 17.2%      |
| 8       | 0.0230064200 | 1.45    | 18.1%      |
| 16      | 0.0204378800 | 1.63    | 10.2%      |

## 7. Наблюдения

- Приватные acc и marker размера \(O(A.\texttt{num\_rows})\) на поток - рост памяти с workers.
- Последовательный merge может ограничивать масштабирование при большом nnz.
- Ускорение \(S > 1\) относительно SEQ наблюдается в pipeline при 8-16 workers
  (до \(S \approx 1.63\)); в task_run при 4-8 workers \(S < 1\).

## 8. Вывод

OMP - статическое разбиение столбцов B, явные shared/default(none), без reduction; слияние CCS
после барьера конца parallel. Наилучший результат на стенде: pipeline, 16 workers (\(S = 1.63\)).
