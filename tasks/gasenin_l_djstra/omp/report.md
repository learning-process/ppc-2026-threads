# Поиск кратчайших путей из одной вершины (алгоритм Дейкстры)

- Студент: Гасенин Леонид Вячеславович, группа 3823Б1ФИ3
- Технология: OMP
- Вариант: 21

## 1. Введение

Данный отчёт описывает реализацию алгоритма Дейкстры с применением технологии OpenMP.
Распараллеливание выполнено на трёх ключевых этапах: поиске вершины с минимальным
расстоянием, релаксации рёбер и финальном суммировании. Последовательная версия (`SEQ`)
используется как эталон для проверки корректности и оценки ускорения.

## 2. Постановка задачи

Дан неориентированный полносвязный граф с $n$ вершинами.

- Вес ребра между вершинами $u$ и $v$ равен $|u - v|$.
- **Входные данные:** количество вершин $n$.
- **Выходные данные:** сумма длин кратчайших путей от вершины 0 до всех остальных.
- **Ограничения:** $n > 0$. Теоретическая сумма: $\frac{n(n-1)}{2}$.

## 3. Схема параллельного алгоритма

Алгоритм работает итерационно. Каждая итерация содержит три параллельные фазы:

1. **Параллельный поиск минимума** (`FindGlobalVertexOMP`): каждый поток независимо
   проходит свою часть вершин (`#pragma omp for nowait`) и сохраняет локальный минимум.
   Барьер (`#pragma omp barrier`) синхронизует записи, затем один поток
   (`#pragma omp single`) сводит локальные минимумы в глобальный.
2. **Параллельная релаксация** (`RelaxEdgesOMP`): обновление расстояний для всех
   непосещённых вершин распараллеливается директивой `#pragma omp parallel for`.
   Запись в `dist[vertex]` безопасна: каждый поток обрабатывает уникальный индекс.
3. **Параллельное суммирование** (`CalculateTotalSumOMP`): финальная сумма расстояний
   вычисляется с редукцией `reduction(+ : total_sum)`.

### 3.1 Ключевой фрагмент — поиск глобальной вершины

Файл: `omp/src/ops_omp.cpp`, функция `FindGlobalVertexOMP`.

```cpp
#pragma omp for nowait
for (int index = 0; index < n; ++index) {
    if (visited[index] == 0 && dist[index] < thread_min) {
        thread_min   = dist[index];
        thread_vertex = index;
    }
}
local_min[thread_id]    = thread_min;
local_vertex[thread_id] = thread_vertex;

#pragma omp barrier

#pragma omp single
{
    InType global_min = inf;
    for (int t = 0; t < num_threads; ++t) {
        if (local_min[t] < global_min) {
            global_min    = local_min[t];
            global_vertex = local_vertex[t];
        }
    }
    if (global_vertex != -1 && global_min != inf)
        visited[global_vertex] = 1;
    else
        global_vertex = -1;
}
```

`nowait` позволяет потокам не ждать конца цикла друг друга; `#pragma omp barrier`
гарантирует завершение всех записей в `local_min` до чтения в `single`.

### 3.2 Ключевой фрагмент — релаксация рёбер

Файл: `omp/src/ops_omp.cpp`, функция `RelaxEdgesOMP`.

```cpp
#pragma omp parallel for default(none) \
    shared(n, inf, global_vertex, dist, visited)
for (int vertex = 0; vertex < n; ++vertex) {
    if (visited[vertex] == 0 && vertex != global_vertex) {
        const InType weight   = std::abs(global_vertex - vertex);
        const InType new_dist = dist[global_vertex] + weight;
        dist[vertex] = std::min(dist[vertex], new_dist);
    }
}
```

Классификация переменных:

- **`shared`:** `dist`, `visited`, `global_vertex`, `n`, `inf`.
- **Приватные (неявно, переменные цикла):** `vertex`, `weight`, `new_dist`.
- **`schedule`:** явно не задан → используется `schedule(static)` по умолчанию.
  OpenMP равномерно делит итерации цикла между потоками блоками одинакового размера.
  Для данной задачи (`n = 200`) статическое разбиение оптимально: нагрузка равномерна.
- **`reduction`:** не используется в `RelaxEdgesOMP` (каждый поток пишет в
  уникальный `dist[vertex]`); редукция применяется в `CalculateTotalSumOMP`:
  `reduction(+ : total_sum)`.

## 4. Детали реализации

- **Заголовочный файл:** `omp/include/ops_omp.hpp`.
- **Файл реализации:** `omp/src/ops_omp.cpp`.
- **Инициализация** (`PreProcessingImpl`): `dist_` и `visited_` выделяются в куче;
  `dist_[0] = 0`, остальные — `INF`. Проверка $n > 0$ — в `ValidationImpl`.
- **Число потоков:** определяется динамически через `GetNumThreads()` (обёртка над
  `omp_get_num_threads()` внутри `#pragma omp single`). Позволяет системе OpenMP
  использовать столько потоков, сколько доступно в данной конфигурации.
- **`default(none)`** во всех параллельных регионах обеспечивает явную декларацию
  переменных и исключает случайный `shared`.
- **Память:** $O(n)$ для `dist_`, $O(n)$ для `visited_`, $O(T)$ для буферов
  локальных минимумов, где $T$ — число потоков. Матрица смежности не создаётся.

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

Корректность OMP-версии проверена функциональными тестами из `tests/functional/main.cpp`
для $n \in \{3, 5, 7\}$. Во всех случаях результат совпал с теоретическим значением
$\frac{n(n-1)}{2}$ и с результатами SEQ-версии.

### 6.2 Производительность

Определения столбцов таблицы:

- **Режим** — сценарий (`task_run` / `pipeline`).
- **MPI-процессов** — число MPI-рангов; в данной технологии это основной параметр
  масштабирования, так как число OMP-потоков фиксировано системой.
- **Время, с** — время выполнения по фреймворку.
- **Ускорение** — $S = T_{\text{seq}} / T_{\text{omp}}$, baseline SEQ `task_run` = 0.01629 с.
- **Эффективность** — не вычисляется, так как число OMP-потоков динамическое (—).

| Режим     | MPI-процессов | Время, с | Ускорение | Эффективность |
|-----------|---------------|----------|-----------|---------------|
| task\_run | 2             | 0.00633  | 2.57×     | —             |
| task\_run | 4             | 0.01987  | 0.82×     | —             |
| task\_run | 6             | 0.01602  | 1.02×     | —             |
| pipeline  | 2             | 0.02262  | 0.72×     | —             |
| pipeline  | 4             | 0.04794  | 0.34×     | —             |
| pipeline  | 6             | 0.07057  | 0.23×     | —             |

**Анализ.** При `task_run` с двумя MPI-процессами OMP демонстрирует ускорение 2.57×.
С ростом числа MPI-процессов производительность падает: 4 и 6 процессов делят одни
и те же 6 физических ядер, каждый запуская собственный пул потоков OpenMP.
Возникает oversubscription — накладные расходы на синхронизацию барьеров
(`#pragma omp barrier`) и переключение контекстов начинают доминировать над
полезной работой. Режим `pipeline` усугубляет картину из-за оверхеда конвейерной
инфраструктуры. Для $n = 200$ этот оверхед особенно заметен.

## 7. Выводы

OMP-реализация корректна и обеспечивает ускорение 2.57× при двух MPI-процессах.
Масштабирование ограничено числом физических ядер и накладными расходами барьеров.
Для максимального ускорения рекомендуется запускать один MPI-процесс и задавать
`OMP_NUM_THREADS` равным числу доступных ядер.

## 8. Список литературы

1. Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
2. OpenMP Architecture Review Board. OpenMP Application Programming Interface. <https://www.openmp.org/specifications/>
3. Документация по курсу. <https://learning-process.github.io/parallel_programming_course/ru/index.html>
