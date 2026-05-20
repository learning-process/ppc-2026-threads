# Поразрядная сортировка для вещественных чисел (тип double) с простым слиянием

- **Студент:** Крымова Кристина Дмитриевна
- **Группа:** 3823Б1Пмоп3
- **Технология:** ALL (MPI + OpenMP)
- **Вариант:** 19

---

## 1. Контекст

Рассматривается гибридная параллельная реализация поразрядной сортировки
массива вещественных чисел типа `double` с использованием:

- распределённой памяти (MPI [1]),
- общей памяти (OpenMP [2]).

**Цель:**

- создание реализации, эффективно использующей технологии `OpenMP` и `MPI`
  совместно,
- масштабирование алгоритма на уровне процессов и потоков.

---

## 2. Постановка задачи

Дан массив `A` размера `n` с элементами типа `double`.
Необходимо отсортировать массив по возрастанию.

### Условия корректности

- Проверяется, что массив не пуст,
- Результат должен быть отсортирован по неубыванию,
- Погрешность: **1e-14**,
- Сравнение с SEQ.

---

## 3. Алгоритм (MPI уровень)

### MPI синхронизация

- `MPI_Bcast` — синхронная рассылка данных (барьерная операция)
- `MPI_Scatterv` — распределение данных с ожиданием завершения
- `MPI_Send/Recv` — блокирующие операции, синхронизирующие процессы
- Неявная синхронизация происходит в момент вызова этих функций

### Общая схема

На уровне процессов алгоритм состоит из нескольких этапов:

1. **Broadcast размера данных** — корневой процесс рассылает всем
   процессам размер массива (`MPI_Bcast`),

2. **Распределение данных (Scatter)** — массив равномерно распределяется
   между процессами с учётом возможного остатка:

```cpp
int chunk = total_size / size_comm;
int rem = total_size % size_comm;
for (int i = 0; i < size_comm; ++i) {
    send_counts[i] = chunk + (i < rem ? 1 : 0);
    offsets[i] = (i == 0) ? 0 : offsets[i - 1] + send_counts[i - 1];
}
```

3. **Локальная сортировка** — каждый процесс сортирует свою часть
   массива с использованием OpenMP для параллельной сортировки блоков,

4. **Сбор результатов (Gather)** — отсортированные части собираются
   на корневом процессе (`MPI_Recv` от каждого процесса),

5. **Слияние** — корневой процесс последовательно сливает
   отсортированные части (простое слияние),

6. **Broadcast результата** — отсортированный массив рассылается
   всем процессам для успешного прохождения тестов.

### Разбиение данных

```cpp
int chunk = total_size / size_comm;
int rem = total_size % size_comm;
send_counts[i] = chunk + (i < rem ? 1 : 0);
```

### Узкие места

- Корневой процесс выполняет:
  - распределение и рассылку данных,
  - последовательное слияние всех частей (O(n) операций),
  - рассылку результата,
- Коммуникации: передача данных между процессами,
- Возможный дисбаланс при неравномерном распределении.

## 4. Алгоритм (OpenMP уровень)

Внутри каждого процесса используется OpenMP для ускорения
локальной сортировки.

### Сортировка блоков

Параллелизм вводится на уровне блоков:

```cpp
#pragma omp parallel for default(none) shared(ull_arr, arr, size)
for (int i = 0; i < size; ++i) {
    ull_arr[i] = DoubleToULL(arr[i]);
}
```

Каждый поток преобразует свою часть массива, затем выполняется
последовательная LSD сортировка (так как для малых блоков overhead
параллелизации не оправдан).

---

## 5. Схема распараллеливания

### MPI уровень

- Процессы независимы,
- Обмен через:
  - `MPI_Bcast` (размер данных и результат),
  - `MPI_Scatterv` (распределение данных),
  - `MPI_Send/Recv` (сбор результатов).

### OpenMP уровень

- Параллелизм внутри каждого процесса при преобразовании
  double → uint64_t,
- `#pragma omp parallel for` для векторизации.

### Разделение данных

- **shared**: `arr`, `ull_arr`,
- **private**: индексы (локальны в цикле).

### Синхронизация

- Отсутствует — каждый поток пишет в свою область памяти,
- Неявный барьер в конце `#pragma omp parallel for`.

---

## 6. Детали реализации

**Файлы:**

- `all/include/ops_all.hpp`
- `all/src/ops_all.cpp`

**Pipeline:**

- **ValidationImpl** — всегда возвращает true,
- **PreProcessingImpl** — подготовка данных,
- **RunImpl**:
  - MPI: `Bcast` размера, `Scatterv` данных,
  - OpenMP: параллельное преобразование double → uint64_t,
  - последовательная LSD сортировка на каждом процессе,
  - `Gather` результатов,
  - слияние на корневом процессе,
  - `Bcast` результата,
- **PostProcessingImpl** — проверка корректности.

### Ключевая идея

- MPI → распределяет данные и собирает результат,
- OpenMP → ускоряет преобразование данных на каждом процессе.

---

## 7. Проверка корректности

Корректность проверялась:

- сравнением с реализацией SEQ на различных данных:
  - граничные случаи,
  - случайные массивы,
  - отсортированные и обратно отсортированные массивы.

Сравнение выполняется с точностью **1e-14**.

Для обеспечения прохождения тестов всеми процессами используется
рассылка результирующего массива всем процессам в `BroadcastResult`.

---

## 8. Экспериментальная среда

- **CPU**: Apple M1 (8 ядер)
- **RAM**: 16 GB
- **OS**: macOS 15.6.1
- **Компилятор**: Clang 21.1.7
- **MPI**: Open MPI v5.0.9
- **OpenMP**: _OPENMP=202011
- **Сборка**: Release

### Сборка и запуск

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON \
  -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cd build/bin

export PPC_NUM_THREADS=<threads>
mpirun -np <num_procs> ./ppc_perf_tests \
  --gtest_filter="*task_run_krymova_k_lsd_sort_merge_double_all_enabled"
```

## 9. Результаты

### Параметры эксперимента

| Параметр | Значение |
| --- | --- |
| Размер массива | 20,000,000 элементов |
| Тип данных | `double` |
| Диапазон значений | случайные числа |

### Результаты измерений

| procs x threads | workers | time (s) | speedup | efficiency |
| --- | --- | --- | --- | --- |
| 1 x 1 | 1 | 0.372 | 1.45 | 1.45 |
| 1 x 2 | 2 | 0.370 | 1.46 | 0.73 |
| 1 x 3 | 3 | 0.371 | 1.45 | 0.48 |
| 1 x 6 | 6 | 0.369 | 1.46 | 0.24 |
| 2 x 1 | 2 | 0.318 | 1.70 | 0.85 |
| 2 x 2 | 4 | 0.319 | 1.69 | 0.42 |
| 2 x 3 | 6 | 0.317 | 1.70 | 0.28 |
| 2 x 6 | 12 | 0.319 | 1.69 | 0.14 |

### Наблюдения

- Использование 2 процессов даёт ускорение **1.70x** — лучше,
  чем 1 процесс
- Добавление потоков внутри процесса не улучшает результат
  из-за memory-bound характера задачи
- Накладные расходы на MPI коммуникации (Scatter/Gather)
  ограничивают масштабируемость
- Лучший результат: **2 процесса × 1 поток = 1.70x ускорения**

---

## 10. Выводы

- Реализован гибридный алгоритм (MPI + OpenMP)
- Достигнуто ускорение до **1.70x**
- Использование двух процессов даёт выигрыш по сравнению с одним
- Дальнейшее увеличение числа потоков не эффективно
  из-за ограничений памяти
- Алгоритм может быть полезен на системах с распределённой памятью,
  где данные не помещаются в память одного узла

---

## Список литературы

1. Документация MPI:
   <https://www.mpich.org>
2. Документация OpenMP:
   <https://www.openmp.org/>
3. LSD Radix Sort:
   <https://en.wikipedia.org/wiki/Radix_sort>

---

## Приложение

### Фрагмент кода — распределение данных и сбор результатов

```cpp
void KrymovaKLsdSortMergeDoubleALL::ComputeDistribution(
    int total_size, int size_comm,
    std::vector<int> &send_counts,
    std::vector<int> &offsets) {

    int chunk = total_size / size_comm;
    int rem = total_size % size_comm;

    for (int i = 0; i < size_comm; ++i) {
        send_counts[i] = chunk + (i < rem ? 1 : 0);
        offsets[i] = (i == 0) ? 0 : offsets[i - 1] + send_counts[i - 1];
    }
}

void KrymovaKLsdSortMergeDoubleALL::ScatterData(
    int rank, const std::vector<int> &send_counts,
    const std::vector<int> &offsets,
    std::vector<double> &local_data) {

    if (rank == 0) {
        MPI_Scatterv(GetInput().data(), send_counts.data(),
                     offsets.data(), MPI_DOUBLE,
                     local_data.data(), send_counts[rank],
                     MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else if (send_counts[rank] > 0) {
        MPI_Scatterv(nullptr, send_counts.data(), offsets.data(),
                     MPI_DOUBLE, local_data.data(),
                     send_counts[rank], MPI_DOUBLE,
                     0, MPI_COMM_WORLD);
    }
}

void KrymovaKLsdSortMergeDoubleALL::GatherResults(
    int rank, int size_comm,
    const std::vector<int> &send_counts,
    std::vector<double> &local_data) {

    if (rank == 0) {
        std::vector<double> result = local_data;
        for (int i = 1; i < size_comm; ++i) {
            if (send_counts[i] > 0) {
                std::vector<double> recv_buf(send_counts[i]);
                MPI_Recv(recv_buf.data(), send_counts[i],
                         MPI_DOUBLE, i, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                result = SimpleMerge(result, recv_buf);
            }
        }
        GetOutput() = std::move(result);
    } else if (send_counts[rank] > 0) {
        MPI_Send(local_data.data(), send_counts[rank],
                 MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}
```
