# Линейная фильтрация блоков – фильтр Гаусса 3×3 — ALL (MPI + threads)

- Student: Москаев Владимир Александрович
- Technology: ALL
- Variant: 26

## 1. Контекст

Гибридная версия сочетает два уровня параллелизма:

- **Межпроцессный**
– с использованием MPI. Каждый MPI-процесс получает подмножество блоков для обработки.
- **Внутрипроцессный**
– внутри каждого MPI-процесса для фильтрации блоков используется std::thread.

Цель – оценить эффективность гибридной схемы.

## 2. Постановка задачи

Полностью совпадает с последовательной версией (см. seq/report.md).
Вход – изображение, выход – отфильтрованное изображение.
Корректность проверяется сравнением с SEQ.

## 3. Базовый алгоритм

Тот же блочный алгоритм, что в SEQ:
разбиение на блоки 64×64, копирование с отступами,
свёртка с ядром Гаусса, запись результата.

## 4. Межпроцессная схема (MPI)

**Роли rank-ов:**

- Все процессы равноправны при вычислениях: каждый процесс обрабатывает свою часть блоков.
- Процесс с рангом 0 выполняет дополнительную роль:
рассылает исходные данные всем процессам и собирает результаты.

**Распределение данных:**

1. Rank 0 рассылает параметры изображения (width, height, channels) через MPI_Bcast.
2. Размер вектора пикселей и сами данные рассылаются через MPI_Bcast.
    Каждый процесс получает полную копию изображения (репликация, а не распределение).
3. Индексы блоков распределяются через MPI_Scatterv.

**Фрагмент кода распределения блоков:**

```cpp
std::vector<int> all(total_blocks);
for (int i = 0; i < total_blocks; ++i) all[i] = i;

std::vector<int> counts(num_procs);
std::vector<int> displs(num_procs);
int off = 0;
for (int proc = 0; proc < num_procs; ++proc) {
    int cnt = per_proc + (proc < rem ? 1 : 0);
    counts[proc] = cnt;
    displs[proc] = off;
    off += cnt;
}
local_blocks.resize(local_cnt);
MPI_Scatterv(all.data(), counts.data(), displs.data(), MPI_INT,
             local_blocks.data(), local_cnt, MPI_INT, 0, MPI_COMM_WORLD);
```

**Сбор результатов:**

- Каждый процесс отправляет свои обработанные данные процессу 0 через MPI_Send.
- Процесс 0 принимает данные от всех процессов через MPI_Recv.
- Затем процесс 0 рассылает собранный результат всем процессам через MPI_Bcast.

**Фрагмент кода сбора результатов:**

```cpp
if (rank == 0) {
    out.resize(total_bytes);
    std::ranges::copy(output, out.begin());
    for (int src = 1; src < num_procs; ++src) {
        if (recv_counts[src] > 0) {
            MPI_Recv(out.data() + displs[src], recv_counts[src], MPI_UNSIGNED_CHAR,
                     src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
} else {
    if (send_count > 0) {
        MPI_Send(output.data(), send_count, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
    }
}
MPI_Bcast(&out_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(out.data(), out_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
```

**Синхронизация:**

- MPI_Bcast – для распространения исходных данных и итогового результата.
- MPI_Scatterv – для распределения индексов блоков.
- MPI_Send/MPI_Recv – для сбора результатов.
- MPI_Barrier не используется,
так как коллективные операции сами обеспечивают синхронизацию.

## 5. Внутрипроцессная схема (std::thread)

Внутри каждого MPI-процесса фильтрация блоков распараллеливается с помощью std::thread.

**Фрагмент кода (ProcessAssignedBlocksParallel):**

```cpp
int num_threads = std::thread::hardware_concurrency();
num_threads = std::min(num_threads, 8);
num_threads = std::min(num_threads, local_cnt);

std::vector<std::vector<uint8_t>> thread_outputs(num_threads);
std::vector<std::thread> threads;

for (int tid = 0; tid < num_threads; ++tid) {
    int blocks_in_thread = blocks_per_thread_base + (tid < blocks_remainder ? 1 : 0);
    int start = (tid * blocks_per_thread_base) + std::min(tid, blocks_remainder);
    
    threads.emplace_back([&, tid, start, blocks_in_thread]() {
        // вычисление размера выходных данных для потока
        // выделение локального буфера
        // обработка блоков
        // сохранение результата в thread_outputs[tid]
    });
}
for (auto &t : threads) t.join();
// объединение результатов из всех потоков
```

**Особенности:**

- Количество потоков ограничено 8 (даже если hardware_concurrency() больше).
- Каждый поток создаёт свой локальный буфер для результатов.
- После завершения всех потоков результаты объединяются в общий буфер.

**Фильтрация одного блока (ProcessOneBlock):**

- Копирование блока с отступами (CopyBlockWithHalo)
- Фильтрация блока (FilterBlock) – внутри также используется std::thread
- Запись результата в выходной буфер

## 6. Детали реализации

**Файлы:** all/include/ops_all.hpp, all/src/ops_all.cpp

**Ключевые функции:**

- BroadcastImageData() – рассылка параметров и данных через MPI_Bcast
- ScatterBlocks() – распределение индексов блоков через MPI_Scatterv
- ProcessAssignedBlocks() – обработка блоков (выбор последовательного или параллельного режима)
- ProcessAssignedBlocksParallel() – параллельная обработка блоков через std::thread
- ProcessOneBlock() – фильтрация одного блока
- FilterBlock() – фильтрация блока с внутренним распараллеливанием по строкам
- GatherAndBroadcastResult() – сбор результатов через MPI_Send/Recv и MPI_Bcast

**Иерархия параллелизма:**

1. MPI-процессы – распределение блоков между процессами
2. std::thread внутри процесса – обработка блоков
3. std::thread внутри FilterBlock – фильтрация строк блока (второй уровень потоков)

## 7. Проверка корректности

Сравнение с SEQ. Функциональные тесты (из tests/functional/main.cpp):

- Тест 1: 2×2 серое, вход [100,150,200,250] → выход [138,163,188,213]
- Тест 2: 3×3 серое, вход [1,2,3,4,5,6,7,8,9] → выход [2,3,4,4,5,6,7,7,8]
- Тест 3: 2×2 RGB, вход 12 чисел → выход 12 чисел
- Тест 4: 1×1 серое, вход [255] → выход [255]

Все тесты пройдены для конфигураций 1×1, 1×4, 2×4, 4×2.

## 8. Экспериментальная среда

- CPU: Intel Core i5-11400H
- RAM: 16 ГБ
- OS: Windows 10
- Компилятор: MSVC 2022 с поддержкой MPI
- Размер задачи: 2048×2048×3 пикселя

**Переменные окружения:** PPC_NUM_THREADS (задаёт число потоков внутри MPI-процесса)

**Команды запуска:**

```bash
$env:PPC_NUM_THREADS=1; mpiexec -n 1 .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_moskaev_v_lin_filt_block_gauss_3_all_enabled"
$env:PPC_NUM_THREADS=4; mpiexec -n 1 .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_moskaev_v_lin_filt_block_gauss_3_all_enabled"
$env:PPC_NUM_THREADS=4; mpiexec -n 2 .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_moskaev_v_lin_filt_block_gauss_3_all_enabled"
$env:PPC_NUM_THREADS=2; mpiexec -n 4 .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_moskaev_v_lin_filt_block_gauss_3_all_enabled"
```

## 9. Результаты

| Ranks | Потоков/rank | Всего workers | Время, с | Speedup | Efficiency |
|-------|--------------|---------------|----------|---------|------------|
| 1     | 1            | 1             | 0.243    | 0.88    | 88%        |
| 1     | 4            | 4             | 0.246    | 0.87    | 22%        |
| 2     | 4            | 8             | 0.165    | 1.30    | 16%        |
| 4     | 2            | 8             | 0.137    | 1.56    | 20%        |

**Комментарий о масштабируемости и узких местах:**

- Лучший результат ALL: 4 процесса × 2 потока = 0.137 с (Speedup 1.56×).
Это значительно хуже, чем у STL (4.86×) и TBB (4.12×) на 8 потоках.

- Конфигурация 1×1 (0.243 с) медленнее SEQ (0.214 с) из-за оверхеда MPI
(даже при одном процессе выполняются MPI_Bcast, MPI_Scatterv, MPI_Send/Recv).

- Конфигурация 1×4 (0.246 с) работает медленнее, чем 1×1,
 что указывает на оверхед от создания потоков внутри MPI-процесса
(PPC_NUM_THREADS=4 создаёт потоки, но при одном процессе это не даёт выигрыша).

- При 2×4 и 4×2 время уменьшается, но speedup всё равно ниже 2×.
Накладные расходы на MPI-коммуникации (рассылка данных, сбор результатов)
 не окупаются на одном узле.

- В коде используется два уровня потоков: ProcessAssignedBlocksParallel создаёт потоки для блоков,
а FilterBlock внутри создаёт потоки для строк.
Это может приводить к избыточному количеству потоков и увеличению оверхеда.

**Примечание по эффективности:**
Efficiency в таблице выше рассчитана как Speedup / (ranks × threads_per_rank),
то есть по общему числу workers. Нормировка выбрана именно по общему числу работников,
чтобы корректно сравнивать ALL с другими backend-ами (OMP, TBB, STL), где workers = число потоков.

## 10. Выводы

ALL версия не даёт выигрыша по сравнению с чистыми STL/TBB на одном компьютере. Основные проблемы:

1. Оверхед на MPI-коммуникации (Bcast, Scatterv, Send/Recv)
2. Репликация данных (каждый процесс получает копию всего изображения)
3. Двухуровневый параллелизм приводит к созданию избыточного количества потоков

MPI имеет смысл использовать только при распределении данных между несколькими узлами,
когда изображение не помещается в память одного узла.
Для одного узла с 8 потоками STL и TBB показывают значительно лучшее ускорение.
