# Умножение плотных матриц. Элементы типа `double`. Алгоритм Штрассена

- Студент: Ахметов Даниил Данисович, группа 3823Б1ПР2
- Технология: SEQ, OMP, TBB, STL, ALL
- Вариант: 3
- Локальные отчёты: `seq/report.md`, `omp/report.md`, `tbb/report.md`,
  `stl/report.md`, `all/report.md`

## 1. Введение

Задача - умножение плотных квадратных матриц типа `double` алгоритмом
Штрассена. Реализованы пять backend-ов: `SEQ`, `OMP`, `TBB`, `STL` и `ALL`.

Цель работы - сравнить разные модели параллелизма на одной и той же
математической задаче. Алгоритм Штрассена удобен для этого, потому что
на верхнем уровне есть семь независимых ветвей `M1..M7`, но при этом
сохраняются накладные расходы на padding, временные буферы и сборку
результата.

Ожидаемый результат: корректное умножение для всех backend-ов
и сравнение производительности относительно последовательного baseline.

## 2. Постановка задачи

Требуется вычислить:

```text
C = A * B
```

где `A` и `B` - квадратные матрицы размера `n × n` с элементами `double`.

**Формат входа** (`InType = std::vector<double>`):

- `input[0]` - размер `n`;
- `input[1 .. n*n]` - матрица `A` (row-major);
- `input[n*n+1 .. 2*n*n]` - матрица `B` (row-major).

**Формат выхода** (`OutType = std::vector<double>`):

- матрица `C` размера `n × n` (row-major).

**Ограничения:**

- `n > 0`;
- `input.size() == 1 + 2 * n * n`;
- если `n` не степень двойки, матрицы дополняются нулями до ближайшей
  степени двойки;
- после вычисления padding удаляется.

**Критерий корректности:**

- результат совпадает с эталонным обычным умножением матриц;
- допуск сравнения: `1e-7`.

## 3. Описание алгоритма (базового/последовательного)

Базовый алгоритм - итеративный Штрассен без рекурсии функций.

**Шаги:**

1. Если `n <= kThreshold` (`kThreshold = 64`), выполняется обычное
   умножение `O(n³)`.
2. Иначе матрицы делятся на квадранты `A11..A22`, `B11..B22`.
3. Вычисляются промежуточные произведения `M1..M7`.
4. Из них собираются блоки `C11..C22`.
5. Квадранты объединяются в итоговую матрицу.

**Практические решения:**

- состояние хранится в стеке структур `Frame`, а не в глубокой рекурсии;
- для малых блоков используется обычное умножение, чтобы снизить
  накладные расходы на временные буферы.

**Асимптотика:**

- время: лучше `O(n³)`, теоретически около `O(n^2.81)`;
- память: выше, чем у наивного умножения, из-за промежуточных буферов.

## 4. Схема распараллеливания

### 4.1 OMP

- **Декомпозиция:** параллелизация внутри итеративного стека, а не только
  на верхнем уровне `M1..M7`.
- **Где используется OpenMP:** `StandardMultiply`, `Split`, `Merge`,
  `AddInto`, `SubInto`, padding и копирование результата.
- **Планирование:** `schedule(static)`, `default(none)`.
- **Условие включения:** `if (size >= kParallelThreshold)`,
  `kParallelThreshold = 256`.
- **Синхронизация:** неявные барьеры в конце `#pragma omp parallel for`;
  `reduction`/`atomic`/`critical` не нужны - потоки пишут в разные элементы.

### 4.2 TBB

- **Декомпозиция:** верхний уровень `M1..M7` - семь независимых задач.
- **Примитив:** `oneapi::tbb::parallel_invoke`.
- **Планирование:** каждая ветвь - отдельная лямбда; `parallel_for`
  не используется, так как задачи дискретные, а не диапазон.
- **Контроль потоков:** `global_control::max_allowed_parallelism`
  по `ppc::util::GetNumThreads()`.
- **Внутренняя часть:** последовательный `StrassenSeqImpl` с блочным
  умножением (`kCutoff = 256`, `kBlockSize = 64`).

### 4.3 STL

- **Декомпозиция:** верхний уровень `M1..M7`.
- **Примитив:** `std::async(std::launch::async, ...)` + `std::future`.
- **Синхронизация:** `future.get()` на главном потоке.
- **Данные:** каждая задача пишет в свой буфер `m1..m7`, гонок нет.

### 4.4 ALL (MPI + потоки)

**MPI (межпроцессный уровень):**

- `MPI_Comm_rank` - определение rank-а;
- каждый rank выполняет полное умножение (без распределения `M1..M7`
  между процессами);
- `MPI_Barrier(MPI_COMM_WORLD)` в конце `RunImpl()`.

**Потоки (внутрипроцессный уровень):**

- основное вычисление - OMP-схема (как в `OMP`);
- после вычисления `C` - демонстрационные блоки:
  - OpenMP (`#pragma omp parallel`, только `rank == 0`);
  - STL (`std::thread` + `atomic`);
  - oneTBB (`tbb::parallel_for`).

Демонстрационные блоки не изменяют результат умножения.

## 5. Детали реализации

**Структура кода:**

| Backend | Заголовок | Исходник |
| --- | --- | --- |
| SEQ | `seq/include/ops_seq.hpp` | `seq/src/ops_seq.cpp` |
| OMP | `omp/include/ops_omp.hpp` | `omp/src/ops_omp.cpp` |
| TBB | `tbb/include/ops_tbb.hpp` | `tbb/src/ops_tbb.cpp` |
| STL | `stl/include/ops_stl.hpp` | `stl/src/ops_stl.cpp` |
| ALL | `all/include/ops_all.hpp` | `all/src/ops_all.cpp` |

**Общий пайплайн задачи** (для всех backend-ов):

- `ValidationImpl()` - проверка формата входа;
- `PreProcessingImpl()` - выделение выходного буфера;
- `RunImpl()` - вычисление;
- `PostProcessingImpl()` - проверка размера результата.

**Крайние случаи:**

- `n <= 64` - сразу обычное умножение;
- `n`, не являющийся степенью двойки - padding и обрезка результата;
- пустой или некорректный вход отклоняется в `ValidationImpl()`.

**Память:**

- матрицы хранятся в `std::vector<double>`;
- в итеративном Штрассене создаются временные буферы для квадрантов
  и `M1..M7` на каждом уровне стека.

## 6. Экспериментальная установка

**Аппаратура и ОС:**

- CPU: AMD Ryzen 7 2700, 8 ядер / 16 потоков;
- RAM: 16 GB;
- OS: Windows 10 x64;
- IDE: Visual Studio Code.

**Инструменты сборки:**

- compiler: MSVC (сборка через CMake);
- build type: `Release`.

**Окружение:**

- `PPC_NUM_THREADS=4` для `SEQ`, `OMP`, `TBB`, `STL` и `ALL` (1 rank);
- `PPC_NUM_THREADS=2`, `mpiexec -n 2` для `ALL`;
- `OMP_NUM_THREADS=4` при запуске OMP/ALL.

**Генерация данных:**

- функциональные тесты: случайные матрицы из `[-10, 10]`,
  размеры `64`, `128`, `256` (`tests/functional/main.cpp`);
- performance-тесты: матрицы `1024 × 1024`, случайные значения
  из `[-10, 10]` (`tests/performance/main.cpp`).

**Команды сборки и запуска:**

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel

cd build\bin
$env:PPC_NUM_THREADS = "4"

# Correctness
.\ppc_func_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*
mpiexec -n 2 .\ppc_func_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*_all_*

# Performance
.\ppc_perf_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*
mpiexec -n 2 .\ppc_perf_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*_all_*
```

**Метрики:**

```text
Speedup = T_seq / T_backend
Efficiency = Speedup / Workers × 100%
```

В таблицу включены значения `task_run`. Режим `pipeline` запускался,
но в сравнение не включался.

## 7. Результаты и обсуждение

### 7.1 Корректность

Корректность проверялась параметризованными функциональными тестами
(`tests/functional/main.cpp`):

- для каждого backend-а генерируются входные матрицы;
- результат сравнивается с наивным эталонным умножением
  в `CheckTestOutputData`;
- допуск: `1e-7`.

Проверенные размеры: `64`, `128`, `256`.

Локально:

- `SEQ`, `OMP`, `TBB`, `STL` - все тесты пройдены;
- `ALL` - пройдены на `mpiexec -n 2`.

### 7.2 Производительность

| Режим | Число потоков | Время, с | Ускорение | Эффективность |
| --- | ---: | ---: | ---: | ---: |
| seq | 1 | 0.938515 | 1.00 | N/A |
| omp | 4 | 1.705478 | 0.55 | 13.8% |
| tbb | 4 | 0.103911 | 9.03 | 225.8% |
| stl | 4 | 0.100045 | 9.38 | 234.5% |
| all (1 rank) | 4 | 1.662727 | 0.56 | 14.1% |
| all (2 ranks × 2 threads) | 4 | 2.187852 | 0.43 | 21.5% |

**Обсуждение:**

- `TBB` и `STL` показали лучшее локальное время. Верхний уровень
  `M1..M7` хорошо ложится на независимые задачи (`parallel_invoke`
  и `std::async`).
- `OMP` и `ALL` оказались медленнее `SEQ` на `1024 × 1024` и 4 потоках.
  Причина - накладные расходы OpenMP на множество относительно небольших
  циклов в итеративной реализации.
- Эффективность выше 100% у `TBB`/`STL` объясняется тем, что это
  другая архитектура верхнего уровня, а не просто добавление потоков
  к той же схеме, что у `SEQ`.
- `ALL` на `mpiexec -n 2` медленнее, потому что каждый rank выполняет
  полное умножение, а не распределяет подзадачи.

Замеры - контрольные локальные прогоны, без серии повторов с медианой.

## 8. Заключение

- Последовательная версия (`SEQ`) работает корректно и пригодна
  как baseline.
- Лучшие результаты по времени на локальных замерах показали `TBB`
  и `STL`.
- `OMP` корректен и используется как основа для `ALL`, но на текущей
  конфигурации не дал ускорения.
- `ALL` выполняет требование семестра потоков: демонстрирует OpenMP,
  STL, TBB и MPI, проходит функциональные тесты на `mpiexec -n 2`.
- Ограничение: без серии повторов и замеров на разных `ranks × threads`
  выводы по масштабируемости предварительные.

## 9. Источники

1. [Parallel Programming Course - threading tasks](https://learning-process.github.io/parallel_programming_course/ru/common_information/threading_tasks.html)
2. V. Strassen, *Gaussian elimination is not optimal*, Numerische Mathematik, 1969.
3. [OpenMP Specifications](https://www.openmp.org/specifications/)
4. [MPI Forum Documentation](https://www.mpi-forum.org/docs/)
5. [oneTBB Documentation](https://uxlfoundation.github.io/oneTBB/)
6. [cppreference: std::async](https://en.cppreference.com/w/cpp/thread/async)

## Приложение (опционально)

Ключевые файлы (пути относительно папки задачи):

- `seq/src/ops_seq.cpp` - baseline, итеративный Штрассен;
- `omp/src/ops_omp.cpp` - OpenMP внутри стека;
- `tbb/src/ops_tbb.cpp` - `parallel_invoke` для `M1..M7`;
- `stl/src/ops_stl.cpp` - `std::async` для `M1..M7`;
- `all/src/ops_all.cpp` - OMP-вычисление + демо технологий + MPI;
- `tests/functional/main.cpp` - проверка корректности;
- `tests/performance/main.cpp` - performance-тесты (`1024 × 1024`).

```cpp
// all/src/ops_all.cpp - фрагмент демонстрации технологий после вычисления C
if (rank == 0) {
  std::atomic<int> counter(0);
#pragma omp parallel default(none) shared(counter) num_threads(ppc::util::GetNumThreads())
  counter++;
}
// ... std::thread и tbb::parallel_for ...
MPI_Barrier(MPI_COMM_WORLD);
```
