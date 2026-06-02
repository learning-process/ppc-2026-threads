# Умножение разреженных матриц (CCS, double) — TBB

- Student: Дилшодов Адхам Умидович, группа 3823Б1ПР4
- Technology: oneTBB (oneAPI Threading Building Blocks)
- Variant: 5

## 1. Контекст

TBB-версия реализует ту же задачно-ориентированную декомпозицию по столбцам результата, что и OMP,
но средствами `tbb::parallel_for` и `tbb::enumerable_thread_specific`.
В отличие от OpenMP, TBB сам управляет work-stealing-ом и временем жизни worker-потоков,
поэтому ручного создания `parallel`-региона не требуется.

## 2. Постановка задачи

Без изменений: `C = A · B`, обе матрицы в CCS, инварианты структуры и численная точность те же,
что в [seq/report.md](../seq/report.md).

## 3. Базовый алгоритм

Идентичен OMP: для каждого столбца `j` матрицы `B` собирается sparse accumulator
(`acc[m]`, `marker[m]`, `used_rows`), затем `used_rows` сортируется и отфильтрованные пары
выгружаются как `j`-й столбец `C`. Финальная склейка локальных результатов в общий CCS —
сериальная фаза `BuildOutputFromColumns`.

## 4. Схема распараллеливания

- **Примитив:** `tbb::parallel_for(blocked_range<int>(0, cols_B), lambda)`. `blocked_range`
  позволяет TBB разбивать диапазон на чанки и переразбивать их при необходимости через
  partitioner по умолчанию (`auto_partitioner`), который комбинирует начальное равномерное
  разбиение со стилом work-stealing.
- **Grainsize:** не указан явно (используется значение по умолчанию). Для текущей задачи
  на 2000 столбцов это разумно: средний чанк — десятки колонок, что покрывает overhead
  на захват задачи и при этом достаточно мелкий для динамической балансировки.
- **Partitioner:** `auto_partitioner` (по умолчанию). Альтернатива — `simple_partitioner`
  с фиксированным grainsize — была бы предпочтительна при сильно перекошенной плотности
  столбцов, но для банд-матриц нагрузка достаточно равномерная.
- **Thread-local буферы:** `tbb::enumerable_thread_specific<ScratchData>` создаёт по одному
  комплекту `acc/marker/used_rows` на worker-поток. Это аналог `private`-буферов в OMP,
  но без явного `parallel`-региона — TBB сам решает, сколько потоков выделить под пул.
- **Конкуренция:** ограничивается рантаймом TBB через
  `tbb::global_control::max_allowed_parallelism`.
  В тестах курса значение берётся из `PPC_NUM_THREADS` через инфраструктуру.

Ключевой фрагмент:

```cpp
// File: tbb/src/ops_tbb.cpp
ScratchData exemplar;
exemplar.acc.assign(lhs.rows_count, 0.0);
exemplar.marker.assign(lhs.rows_count, -1);

tbb::enumerable_thread_specific<ScratchData> tls(exemplar);

tbb::parallel_for(
    tbb::blocked_range<int>(0, rhs.cols_count),
    [&](const tbb::blocked_range<int> &range) {
      auto &scratch = tls.local();
      for (int rhs_col = range.begin(); rhs_col < range.end(); ++rhs_col) {
        AccumulateColumnProduct(lhs, rhs, rhs_col, scratch, column_results[rhs_col]);
      }
    });
```

Расшифровка:

- `exemplar` хранит «прообраз» thread-local буфера.
  `enumerable_thread_specific` лениво клонирует его при первом обращении на каждом потоке;
- `tls.local()` возвращает ссылку на буфер текущего потока за O(1) (по сути TLS lookup);
- `blocked_range` сам решает, как делить диапазон, и TBB при свободном потоке может
  «украсть» половину чанка у занятого.

## 5. Детали реализации

**Файлы:** [tbb/include/ops_tbb.hpp](include/ops_tbb.hpp), [tbb/src/ops_tbb.cpp](src/ops_tbb.cpp).

Локальный struct `ScratchData` группирует `acc`, `marker`, `used_rows` в один объект,
чтобы передавать по ссылке без переинициализации. Между итерациями `acc` не обнуляется —
это безопасно, потому что `marker[i] == j` обновляется при первом обращении в новой колонке,
и старое значение `acc[i]` затирается без чтения.

Объединение результата (`BuildOutputFromColumns`) — последовательное, по тому же шаблону
что и в OMP: префиксная сумма `col_ptrs` и последовательная конкатенация.
Параллелить эту фазу через `parallel_scan` нецелесообразно: на 2000 столбцов время сборки
исчезающе мало по сравнению с основным вычислением.

## 6. Проверка корректности

Все три функциональных теста (`TwoByTwoBasic`, `ThreeByThreeSparse`, `RectangularCheck`)
проходят, выход побитно совпадает с эталоном. Производительный тест на матрице 2000×2000
проходит для всех значений `PPC_NUM_THREADS=1,2,4,8` в обоих режимах (`pipeline`, `task_run`)
— это подтверждает отсутствие гонок и стабильную корректность под разной степенью конкуренции.

Структурно гонок нет: `column_results[j]` адресуется уникальным индексом, `tls.local()`
возвращает приватный буфер, чтения `lhs`/`rhs` идут только в режиме read.

## 7. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | AMD Ryzen 5 5600 (6 ядер / 12 потоков SMT) |
| RAM | 15 GiB |
| OS | Ubuntu 24.04.2 LTS в Docker (WSL2 ядро 6.6.114.1, хост Windows 11) |
| Compiler | g++ Ubuntu 13.x, `-O3` |
| oneTBB | поставляется через `ppc_onetbb` из репозитория курса |
| `PPC_NUM_THREADS` | 1, 2, 4, 8 |
| Размер задачи | A, B — банд 2000×2000, ширина ±50 |

```bash
for N in 1 2 4 8; do
  PPC_NUM_THREADS=$N OMP_NUM_THREADS=$N \
    ./build_dilshodov/bin/ppc_perf_tests \
    --gtest_filter="*Dilshodov*tbb*"
done
```

## 8. Результаты

| Workers | TBB, task_run | TBB, pipeline | Speedup vs SEQ@1 | Speedup vs TBB@1 | Efficiency vs TBB@1 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.0288 | 0.0336 | 6.78 | 1.00 | 100% |
| 2 | 0.0164 | 0.0220 | 11.91 | 1.76 | 88% |
| 4 | 0.0097 | 0.0138 | 20.14 | 2.97 | 74% |
| 8 | 0.0068 | 0.0129 | 28.74 | 4.24 | 53% |

**Комментарий.**

- TBB@1 чуть медленнее OMP@1 (0.0288 vs 0.0246 с) — это типичный overhead на инициализацию
  task-arena и `enumerable_thread_specific`.
- На 4 потоках TBB и OMP сравниваются по времени; на 8 потоках **TBB обходит OMP**
  (0.0068 vs 0.0077 с). Видимо, work-stealing при SMT-нагрузке справляется лучше
  с неравномерностью, чем `schedule(dynamic)` OMP.
- Эффективность относительно TBB@1: 88%/74%/53% — заметно лучше, чем у OMP (80%/69%/40%).
  Это согласуется с тем, что TBB разработан под задачно-ориентированный паттерн
  и эффективнее адаптируется к фоновой нагрузке системы.

## 9. Выводы

TBB показывает чуть более высокий overhead на одном потоке, но лучше масштабируется
при росте конкуренции, особенно на SMT. Для текущей задачи
`parallel_for + enumerable_thread_specific + blocked_range` без явного grainsize
и без `simple_partitioner` оказались оптимальной комбинацией: код короткий,
нет ручного управления потоками, накладные расходы амортизируются на 2000 итерациях.

Возможные направления улучшения: вынести `BuildOutputFromColumns` в `parallel_for`
со вторым проходом (после префиксной суммы) — но на текущем размере задачи
это даст менее 5% выигрыша.
