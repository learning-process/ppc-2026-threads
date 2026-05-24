# TBB: Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – столбцовый (CCS).

Обзор задачи, метрики и сборка: [../report.md](../report.md).

- Студент: Борунов Владислав Алексеевич
- Группа: 3823Б1ПР3
- Вариант: 7
- Задача: Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – столбцовый (CCS).
- Технологии: TBB

## 1. Назначение

BorunovVComplexCcsTBB - параллельное умножение через Intel oneTBB: та же логика столбца, что
у OMP (ProcessColumn), оркестрация через task_arena и parallel_for.

## 2. Распараллеливание

Декомпозиция по индексу worker tid - диапазон столбцов j [jstart, jend) как в OMP.

Фаза 1 - вычисление (tbb/src/ops_tbb.cpp):

```cpp
tbb::task_arena arena(num_threads);
arena.execute([&] {
  tbb::parallel_for(
      tbb::blocked_range<int>(0, num_threads, 1),
      [&](const tbb::blocked_range<int> &r) { ... },
      tbb::static_partitioner());
});
```

Фаза 2 — копирование в C.values / C.row_indices (после расчёта col_ptrs и thread_offsets на главном потоке):

```cpp
tbb::parallel_for(
    tbb::blocked_range<int>(0, num_threads, 1),
    [&](const tbb::blocked_range<int> &r) { std::copy(...); },
    tbb::static_partitioner());
```

## 3. Синхронизация, blocked_range, partitioner

| Параметр      | Значение             | Смысл                                                      |
|---------------|----------------------|------------------------------------------------------------|
| blocked_range | [0, num_threads)     | один чанк = один tid                                       |
| grainsize     | 1                    | минимальный поддиапазон, задача на каждый tid              |
| partitioner   | static_partitioner() | фиксированное разбиение, предсказуемая привязка к столбцам |

| Механизм                          | Роль                                    |
|-----------------------------------|-----------------------------------------|
| task_arena(num_threads)           | число worker TBB = PPC_NUM_THREADS      |
| t_values[tid], t_row_indices[tid] | нет записи в общий C в фазе 1           |
| acc, marker, touched              | локальные в задаче (аналог private OMP) |
| merge col_ptrs                    | главный поток между execute             |
| вторая parallel_for               | копирование в непересекающиеся сегменты |

mutex не используется.

## 4. Отличия от других реализаций

| Аспект           | TBB                          | OMP                     | SEQ / STL             |
|------------------|------------------------------|-------------------------|-----------------------|
| Оркестрация      | task_arena + 2x parallel_for | один parallel           | - / std::thread       |
| Merge значений   | параллельный copy            | последовательный insert | - / insert после join |
| Алгоритм столбца | Shell sort, marker           | то же                   | sort + bool (SEQ/STL) |

## 5. Конфигурация workers

workers = ppc::util::GetNumThreads() (PPC_NUM_THREADS).

## 6. Результаты (локальный прогон)

Baseline: SEQ task_run = 0.0408020800 с., SEQ pipeline = 0.0332925600 с.

### 6.1 task_run

| Workers | time, с      | speedup | efficiency |
|---------|--------------|---------|------------|
| 4       | 0.0389541200 | 1.05    | 26.2%      |
| 8       | 0.0473973000 | 0.86    | 10.8%      |
| 16      | 0.0371147000 | 1.10    | 6.9%       |

### 6.2 pipeline

| Workers | time, с      | speedup | efficiency |
|---------|--------------|---------|------------|
| 4       | 0.0512741400 | 0.65    | 16.2%      |
| 8       | 0.0334311600 | 1.00    | 12.4%      |
| 16      | 0.0525265200 | 0.63    | 4.0%       |

## 7. Наблюдения

- Лучшие \(S \ge 1\) в режиме task_run (4 и 16 workers).
- В pipeline при 16 workers время растёт (\(S = 0.63\)) - возможное влияние конкуренции за ресурсы с другими этапами конвейера.
- static_partitioner оправдан: нагрузка по tid отличается не более чем на один столбец.

## 8. Вывод

TBB использует task_arena, blocked_range с grainsize=1 и static_partitioner; конкуренция снята
раздельными буферами и двухфазным слиянием. Наилучший результат на стенде: task_run, 16
workers (\(S = 1.10\)) и task_run, 4 workers (\(E = 26.2\%\)).
