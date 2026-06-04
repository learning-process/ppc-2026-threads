# Поиск кратчайших путей из одной вершины (алгоритм Дейкстры)

- Студент: Гасенин Леонид Вячеславович, группа 3823Б1ФИ3
- Технология: ALL (MPI + OpenMP)
- Вариант: 21

## 1. Введение

Данный отчёт описывает гибридную реализацию алгоритма Дейкстры, совмещающую MPI
для межпроцессного параллелизма и OpenMP для внутрипроцессного.
MPI распределяет вершины графа между рангами, каждый из которых независимо обрабатывает
свой блок; OpenMP распараллеливает поиск локального минимума и релаксацию внутри
каждого ранга. Последовательная версия (`SEQ`) служит эталоном корректности и ускорения.

## 2. Постановка задачи

Дан неориентированный полносвязный граф с $n$ вершинами.

- Вес ребра между вершинами $u$ и $v$ равен $|u - v|$.
- **Входные данные:** количество вершин $n$.
- **Выходные данные:** сумма длин кратчайших путей от вершины 0 до всех остальных.
- **Ограничения:** $n > 0$. Теоретическая сумма: $\frac{n(n-1)}{2}$.

## 3. Схема параллельного алгоритма

### 3.1 Конфигурация рангов и потоков

Реализация использует двухуровневый параллелизм: **MPI-ранги × OpenMP-потоки**.

- **MPI-уровень:** $P$ рангов, каждый хранит фрагмент `dist_` и `visited_`
  размером $\approx n/P$.
- **OpenMP-уровень:** внутри каждого ранга запускается $T$ потоков, определяемых
  динамически через `omp_get_num_threads()`. На машине с 6 логическими ядрами при
  запуске `mpiexec -n 2` типичная конфигурация: **2 ранга × 3 потока** = 6 рабочих
  единиц. При `-n 6`: **6 рангов × 1 поток** (oversubscription).

### 3.2 Распределение данных (MPI-уровень)

В `PreProcessingImpl` вершины делятся между `size_` рангами блочным разбиением с
остатком:

- Ранги с номером `< rem` получают блок размером `chunk + 1`.
- Остальные ранги получают блок размером `chunk`.

Каждый ранг хранит только свой фрагмент. Вершина-источник (0) инициализируется
нулём тем рангом, в чей диапазон она попадает.

### 3.3 Итерация алгоритма (два уровня параллелизма)

Каждая из $n$ итераций включает четыре шага:

1. **OpenMP-поиск локального минимума** (`FindThreadLocalMinima`): потоки OpenMP
   проходят свои части локального блока (`#pragma omp for nowait`) и записывают
   результаты в `thread_mins[tid]` / `thread_vertices[tid]`.
2. **Свёртка потоковых минимумов** (`ReduceThreadMinima`): последовательная свёртка
   локальных результатов потоков в пару `(local_min, local_vertex)`.
3. **MPI-редукция глобального минимума** (`MPI_Allreduce` с `min_pair_op`): все ранги
   согласовывают глобальную вершину. `MPI_Allreduce` — коллективная блокирующая
   операция: она неявно синхронизирует все ранги, то есть ни один ранг не продолжает
   итерацию до получения глобального результата.
4. **OpenMP-релаксация** (`UpdateLocalDistances`): `#pragma omp parallel for` обновляет
   расстояния до непосещённых локальных вершин.

По окончании всех итераций: каждый ранг считает локальную сумму (OpenMP `reduction`),
затем `MPI_Allreduce(MPI_SUM)` собирает глобальный результат.

### 3.4 Ключевой фрагмент — пользовательская операция MPI

Файл: `all/src/ops_all.cpp`, функция `MinPairImpl`.

```cpp
void MinPairImpl(void *in, void *inout, const int *len, MPI_Datatype * /*dtype*/) {
  auto *a = static_cast<InType *>(in);
  auto *b = static_cast<InType *>(inout);
  for (int i = 0; i < *len; i += 2) {
    if (a[i] < b[i]) {
      b[i]     = a[i];
      b[i + 1] = a[i + 1];
    }
  }
}
```

Буфер трактуется как массив пар `(расстояние, индекс_вершины)`. При `MPI_Allreduce`
с `count = 1` все ранги получают одну глобальную пару-победитель — вершину с
наименьшим расстоянием.

### 3.5 Ключевой фрагмент — OpenMP-поиск локального минимума

Файл: `all/src/ops_all.cpp`, функция `FindThreadLocalMinima`.

```cpp
#pragma omp parallel default(none) \
    shared(local_n, start_v, dist, visited, thread_mins, thread_vertices, inf)
{
  const int tid = omp_get_thread_num();
  InType t_min = inf;
  InType t_v   = -1;

#pragma omp for nowait
  for (int i = 0; i < local_n; ++i) {
    if (visited[i] == 0 && dist[i] < t_min) {
      t_min = dist[i];
      t_v   = start_v + i;
    }
  }

  thread_mins[tid]     = t_min;
  thread_vertices[tid] = t_v;
}
```

**Классификация переменных:**

- **`shared`:** `local_n`, `start_v`, `dist`, `visited`, `thread_mins`,
  `thread_vertices`, `inf`.
- **Приватные (неявно):** `tid`, `t_min`, `t_v`, `i`.
- **`schedule`:** не задан → `schedule(static)` по умолчанию.
- Запись в `thread_mins[tid]` безопасна: каждый поток использует уникальный индекс.

### 3.6 Ключевой фрагмент — OpenMP-релаксация рёбер

Файл: `all/src/ops_all.cpp`, функция `UpdateLocalDistances`.

```cpp
#pragma omp parallel for default(none) \
    shared(local_n, start_v, dist, visited, global_vertex, global_min)
for (int i = 0; i < local_n; ++i) {
  if (visited[i] == 0) {
    const InType global_i  = start_v + i;
    if (global_i != global_vertex) {
      const InType weight   = std::abs(global_vertex - global_i);
      const InType new_dist = global_min + weight;
      dist[i] = std::min(new_dist, dist[i]);
    }
  }
}
```

**Классификация переменных:**

- **`shared`:** `local_n`, `start_v`, `dist`, `visited`, `global_vertex`, `global_min`.
- **Приватные (неявно):** `i`, `global_i`, `weight`, `new_dist`.
- Каждый поток обрабатывает уникальный индекс `i` — гонок при записи в `dist[i]` нет.

## 4. Детали реализации

- **Заголовочный файл:** `all/include/ops_all.hpp`.
- **Файл реализации:** `all/src/ops_all.cpp`.
- **Инициализация:** `PreProcessingImpl` — блочное разбиение; `dist_[0 - start_v_] = 0`
  на ранге-владельце вершины 0, остальные — `INF`.
- **Число потоков OpenMP:** определяется динамически в начале `RunImpl`:
  `#pragma omp single { num_threads = omp_get_num_threads(); }`.
- **Пользовательская операция MPI:** `MPI_Op_create` в каждом `RunImpl`,
  `MPI_Op_free` до выхода — утечек нет.
- **MPI-синхронизация:** `MPI_Allreduce` — коллективная блокирующая операция,
  неявно выступает барьером для всех рангов на каждой итерации.
- **OpenMP-синхронизация:** неявный барьер в конце `#pragma omp parallel` гарантирует
  завершение `FindThreadLocalMinima` до передачи результата в `MPI_Allreduce`.
- **Память:** $O(n / P)$ для `dist_` и `visited_` на каждом ранге; $O(T)$ для
  буферов потоков. Матрица смежности не создаётся.

## 5. Экспериментальная установка

- **Hardware/OS:**
  - CPU: Intel Core i5-8400 2.80 GHz (6 ядер, 6 потоков)
  - RAM: 8 ГБ
  - OS: Windows 10
- **Toolchain:**
  - IDE: Visual Studio Code
  - Компилятор: GCC
  - Система сборки: CMake
  - Система контроля версий: Git
- **Размер задачи:** $n = 200$ вершин.
- **Команды сборки:**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

- **Команда запуска функциональных тестов:**

```bash
mpiexec -n 2 .\build\bin\ppc_func_tests.exe
```

- **Команды запуска тестов производительности:**

```bash
mpiexec -n 2 .\build\bin\ppc_perf_tests.exe --gtest_filter="Djkstra*"
mpiexec -n 4 .\build\bin\ppc_perf_tests.exe --gtest_filter="Djkstra*"
mpiexec -n 6 .\build\bin\ppc_perf_tests.exe --gtest_filter="Djkstra*"
```

- **Сценарии замера:** `pipeline` и `task_run`.
- **Базовая линия SEQ** (`task_run`): `0.01629 с`.

## 6. Результаты

### 6.1 Корректность

ALL-версия прошла все функциональные тесты из `tests/functional/main.cpp` для
$n \in \{3, 5, 7\}$. Результаты совпали с теоретическим значением $\frac{n(n-1)}{2}$
и с SEQ-версией. Тесты запускались при 2, 4 и 6 MPI-процессах.

### 6.2 Производительность

Определения столбцов таблицы:

- **Режим** — сценарий (`task_run` / `pipeline`).
- **MPI-процессов** — число рангов $P$; фактическое число рабочих единиц ≈ $P \times T$,
  где $T$ — число OMP-потоков на ранг (определяется системой).
- **Время, с** — время выполнения по фреймворку.
- **Ускорение** — $S = T_{\text{seq}} / T_{\text{all}}$, baseline SEQ = 0.01629 с.
- **Эффективность** — не вычисляется ($P \times T$ динамически; —).

| Режим     | MPI-процессов | Время, с | Ускорение |
|-----------|---------------|----------|-----------|
| task\_run | 2             | 0.00448  | 3.64×     |
| task\_run | 4             | 0.08546  | 0.19×     |
| task\_run | 6             | 0.26532  | 0.06×     |
| pipeline  | 2             | 0.02010  | 0.81×     |
| pipeline  | 4             | 0.34253  | 0.05×     |
| pipeline  | 6             | 1.03576  | 0.02×     |

**Анализ.** При `task_run` с двумя MPI-процессами (конфигурация ~2×3 = 6 рабочих
единиц) гибридная версия даёт 3.64× — наивысший результат среди всех технологий
в данном сценарии. С ростом числа рангов до 4 и 6 производительность резко падает:
`MPI_Allreduce` с пользовательской операцией создаёт коммуникационный оверхед,
превышающий вычислительную работу при $n = 200$. При 6 рангах возникает
oversubscription (6 рангов × ≥1 OMP-поток > 6 ядер), что усугубляет ситуацию.
Режим `pipeline` добавляет оверхед конвейерной инфраструктуры.

## 7. Выводы

Гибридная ALL-реализация корректна и при конфигурации **2 ранга × T OMP-потоков**
обеспечивает наивысшее ускорение среди всех технологий (3.64×). При большем числе
рангов коммуникационный оверхед `MPI_Allreduce` перевешивает вычислительный
выигрыш для $n = 200$. Оптимальная конфигурация: один ранг на NUMA-узел с
`OMP_NUM_THREADS` равным числу ядер на узел.

## 8. Список литературы

1. Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
2. OpenMP Architecture Review Board. OpenMP Application Programming Interface. <https://www.openmp.org/specifications/>
3. MPI Forum. MPI: A Message-Passing Interface Standard. <https://www.mpi-forum.org/docs/>
4. Документация по курсу. <https://learning-process.github.io/parallel_programming_course/ru/index.html>
