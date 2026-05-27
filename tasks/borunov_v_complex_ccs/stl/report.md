# STL: Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – столбцовый (CCS)

Обзор задачи, метрики и сборка: [../report.md](../report.md).

- Студент: Борунов Владислав Алексеевич
- Группа: 3823Б1ПР3
- Вариант: 7
- Задача: Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – столбцовый (CCS).
- Технологии: STL

## 1. Назначение

BorunovVComplexCcsSTL - параллельное умножение на std::thread:
разбиение столбцов B, локальные буферы, join() до merge в главном потоке.

## 2. Распараллеливание

разбиение столбцов (как в OMP):

```cpp
start_col = (num_cols * thread_id) / num_threads
end_col   = (num_cols * (thread_id + 1)) / num_threads
```

Число потоков (stl/src/ops_stl.cpp):

```cpp
unsigned int num_threads = std::thread::hardware_concurrency();
if (num_threads == 0) num_threads = 4;
if (num_threads > num_cols) num_threads = num_cols;
```

WorkerThread: для каждого j в диапазоне - накопление (как SEQ),
std::ranges::sort, запись в thread_val, thread_row_idx, thread_col_ptr; порог |z| > 1e-9.

## 3. Синхронизация

```cpp
for (unsigned int i = 0; i < num_threads; ++i) {
  threads.emplace_back(worker, i);
}
for (auto &t : threads) {
  t.join();
}
// merge - только после join
```

Главный поток не пишет в C во время работы workers -
data race отсутствуют. Merge: сборка col_ptrs, insert в C.values / C.row_indices.

## 4. Отличия от других реализаций

| Аспект      | STL                         | OMP / TBB       | SEQ               |
|-------------|-----------------------------|-----------------|-------------------|
| workers     | hardware_concurrency()      | PPC_NUM_THREADS | 1                 |
| Маркер      | `vector<bool>`              | marker[i]!=j    | `vector<bool>`    |
| Сортировка  | std::ranges::sort           | Shell sort      | std::ranges::sort |
| Пул потоков | создаётся на каждый RunImpl | runtime OMP/TBB | —                 |

## 5. Конфигурация workers

workers = std::thread::hardware_concurrency() (не PPC_NUM_THREADS).

На стенде замеров при PPC_NUM_THREADS = 4, 8, 16

## 6. Результаты (локальный прогон)

Baseline: SEQ task_run = 0.0408020800 с., SEQ pipeline = 0.0332925600 с.

### 6.1 task_run

| Workers | time, с      | speedup | efficiency |
|---------|--------------|---------|------------|
| 4       | 0.0579572800 | 0.70    | 17.6%      |
| 8       | 0.0591825200 | 0.69    | 8.6%       |
| 16      | 0.0493978000 | 0.83    | 5.2%       |

### 6.2 pipeline

| Workers | time, с      | speedup | efficiency |
|---------|--------------|---------|------------|
| 4       | 0.0352535600 | 0.94    | 23.6%      |
| 8       | 0.0302626200 | 1.10    | 13.8%      |
| 16      | 0.0333756400 | 1.00    | 6.2%       |

## 7. Наблюдения

- Создание/уничтожение потоков на каждый RunImpl добавляет накладные расходы в task_run.
- pipeline при 8 workers даёт \(S = 1.10\) - сопоставимо с TBB при том же режиме.
- При сравнении на другой машине нужно явно фиксировать
  hardware_concurrency(), иначе efficiency по PPC_NUM_THREADS будет некорректна.

## 8. Вывод

STL - явное разбиение столбцов B, параллельная фаза с обязательным
join() до merge. Наилучший результат на стенде: pipeline, 8 workers (\(S = 1.10\)).
