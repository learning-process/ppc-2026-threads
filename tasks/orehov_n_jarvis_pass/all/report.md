# <Построение выпуклой оболочки – проход Джарвиса> — ALL
- Student: Орехов Никита Антонович
- Technology: ALL
- Variant: 23

## 1. Контекст
Гибридная версия сочетает два уровня параллелизма:
- **Межпроцессный** – с использованием MPI (Message Passing Interface). Каждый MPI-процесс получает подмножество точек, выполняет локальный поиск следующей вершины оболочки, затем результаты собираются на нулевом процессе.
- **Внутрипроцессный** – внутри каждого MPI-процесса для обработки локального диапазона точек используется OpenMP, что позволяет эффективно задействовать все ядра одного узла.

Цель – оценить, даёт ли гибридная схема дополнительный выигрыш по сравнению с чистыми реализациями.

## 2. Постановка задачи
Полностью совпадает с последовательной версией (см. `seq/report.md`). Вход – множество точек, выход – вершины выпуклой оболочки в порядке обхода против часовой стрелки. Обработка крайних случаев (1, 2 точки, дубликаты, коллинеарность) аналогична SEQ.

## 3. Базовый алгоритм
Строим выпуклую оболочку множества точек, последовательно находя следующую вершину с максимальным левым поворотом относительно текущей, начиная с самой левой нижней точки.
Подробности – в отчёте SEQ.

## 4. Межпроцессная схема
**Роли rank-ов**:
- Все процессы равноправны при вычислениях: каждый процесс обрабатывает свою часть массива точек.
- Процесс с рангом 0 выполняет дополнительную роль: собирает глобальные результаты, удаляет дубликаты на этапе предобработки, рассылает данные всем процессам.

**Распределение данных**:
1. Процесс 0 удаляет дубликаты через `std::set`, формирует итоговый вектор `input_`.
2. Размер вектора `vec_size` рассылается всем процессам через `MPI_Bcast`.
3. Сам вектор (координаты точек) упаковывается в плоский массив `double` и также рассылается через `MPI_Bcast`. Каждый процесс восстанавливает свою копию `input_`.

**Параллельный поиск следующей вершины** (функция `FindNext`):
- Каждый процесс получает свою непрерывную часть массива `input_` на основе линейного разбиения:  
  `chunk = n / size`, `rest = n % size`.  
  Для rank `r`: `start = r*chunk + min(r, rest)`, `end = start + chunk + (r < rest ? 1 : 0)`.
- Процесс выполняет локальный поиск лучшей точки в своём диапазоне, используя OpenMP (см. следующий раздел).
- Локальные лучшие точки собираются на процессе 0 через `MPI_Gather`. Каждый процесс отправляет 3 числа: `(x, y, valid_flag)`.
- Процесс 0 выполняет глобальную редукцию (`GlobalReduce`), выбирая лучшую точку среди всех.
- Результат (координаты и флаг) рассылается всем процессам через `MPI_Bcast`.

**Синхронизация**:
- `MPI_Bcast` – для распространения исходных данных и итоговой лучшей точки.
- `MPI_Gather` – для сбора локальных лучших точек на нулевой процесс.
- Явный барьер не используется, так как `MPI_Gather` и `MPI_Bcast` сами обеспечивают синхронизацию.

## 5. Внутрипроцессная схема
Внутри каждого MPI-процесса поиск лучшей точки в выделенном диапазоне распараллеливается с помощью OpenMP. Это реализовано в функции `LocalFindBest`:

- Создаётся параллельная область `#pragma omp parallel default(none) shared(...)`.
- Каждый поток получает свой индекс `tid` и локальную структуру `thread_best`.
- Цикл по элементам диапазона распределяется между потоками с помощью `#pragma omp for nowait`.
- После завершения цикла каждый поток обновляет свою локальную лучшую точку.
- Затем с помощью `#pragma omp critical` потоки объединяют свои результаты в общую `local_best` (поскольку запись в одну переменную требует взаимного исключения).


## 6. Детали реализации
'''cpp
OrehovNJarvisPassALL::BestState OrehovNJarvisPassALL::LocalFindBest(const Point &current, size_t start,
                                                                    size_t end) const {
  BestState local_best;
  // Сохраняем this в константную переменную для shared
  const auto *self = this;
#pragma omp parallel default(none) shared(local_best, start, end, current, self)
  {
    BestState thread_best;
#pragma omp for nowait
    for (size_t i = start; i < end; ++i) {
      const Point &p = self->input_[i];
      if (p == current) {
        continue;
      }
      if (!thread_best.valid || IsBetterPoint(current, p, thread_best.point)) {
        thread_best.point = p;
        thread_best.valid = true;
      }
    }
#pragma omp critical
    {
      if (thread_best.valid) {
        if (!local_best.valid) {
          local_best = thread_best;
        } else {
          local_best = ReduceBestStates(local_best, thread_best, current);
        }
      }
    }
  }
  return local_best;
}


Point OrehovNJarvisPassALL::FindNext(Point current) const {
  const size_t n = input_.size();
  if (n == 0) {
    return current;
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t chunk = n / size;
  size_t rest = n % size;
  size_t start = (static_cast<size_t>(rank) * chunk) + std::min(static_cast<size_t>(rank), rest);
  size_t end = start + chunk + (std::cmp_less(static_cast<size_t>(rank), rest) ? 1 : 0);

  BestState local_best = LocalFindBest(current, start, end);

  std::array<double, 3> local_data = {local_best.valid ? local_best.point.x : 0.0,
                                      local_best.valid ? local_best.point.y : 0.0, local_best.valid ? 1.0 : 0.0};

  std::vector<double> all_data(static_cast<size_t>(size) * 3);
  MPI_Gather(local_data.data(), 3, MPI_DOUBLE, all_data.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::array<double, 3> global_data = {0.0, 0.0, 0.0};
  if (rank == 0) {
    BestState global_best = GlobalReduce(all_data, size, current);
    if (global_best.valid) {
      global_data[0] = global_best.point.x;
      global_data[1] = global_best.point.y;
      global_data[2] = 1.0;
    }
  }

  MPI_Bcast(global_data.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  BestState final_best = FinalizeBestPoint(global_data.data());
  return final_best.valid ? final_best.point : current;
}
'''

## 7. Проверка корректности
- **Сравнение результатов**: выходной вектор `res_` сравнивается с эталонным вектором, заданным в функциональных тестах. Порядок точек важен: оболочка должна быть перечислена, начиная с самой левой нижней точки, против часовой стрелки.

- **Функциональные тесты**:
  - Тест 1: 7 точек, образующих выпуклый многоугольник. Ожидается 5 вершин: (0,0), (0,3), (2,4), (4,2), (3,0).
  - Тест 2: одна точка (0,0) → оболочка {(0,0)}.
  - Тест 3: две точки (0,0) и (1,0) → оболочка {(0,0), (1,0)}.

## 8. Экспериментальная среда
- **CPU**: Intel(R) Core(TM) i5-9300H CPU @ 2.40GHz
- **RAM**: 16,0 ГБ
- **OS**: Windows 10
- **Компилятор**: Visual Studio Community 2022 Release - amd64
- **CMake build type**: `Release`
- **Команда запуска**:

# 1 процесс, 1 поток
set PPC_NUM_THREADS=1 && mpiexec -n 1 ./ppc_perf_tests.exe --gtest_filter="*all*"
# 1 процесс, 4 потока
set PPC_NUM_THREADS=4 && mpiexec -n 1 ./ppc_perf_tests.exe --gtest_filter="*all*"
# 2 процесса, 1 поток на процесс
set PPC_NUM_THREADS=1 && mpiexec -n 2 ./ppc_perf_tests.exe --gtest_filter="*all*"
# 2 процесса, 4 потока на процесс
set PPC_NUM_THREADS=4 && mpiexec -n 2 ./ppc_perf_tests.exe --gtest_filter="*all*"
# 4 процесса, 1 поток на процесс
set PPC_NUM_THREADS=1 && mpiexec -n 4 ./ppc_perf_tests.exe --gtest_filter="*all*"

## 9. Результаты
Базовое время последовательной версии T_seq = 0.18437 с (из отчёта SEQ). Для ALL-версии использовался режим task_run. 

MPI-процессов |	Потоков на процесс | Общее число работников	| Время, с |Speedup|Efficiency|
        1     |         1          |          1             |  0.09590 |  1.92 |   1.92   |
        1     |         4          |          4             |  0.09290 |  1.98 |   0.50   |
        2     |         1          |          2             |  0.20710 |  0.89 |   0.45   |
        2     |         4          |          8             |  0.20682 |  0.89 |   0.11   |
        4     |         1          |          4             |  0.45035 |  0.41 |   0.10   |


## 10. Выводы
Когда гибридная схема оправдана:

 - Режим одного MPI-процесса с OpenMP (фактически OpenMP-версия). В этом случае достигается наилучшее время выполнения: 0.0929 с при 1 процессе и 4 потоках, что соответствует ускорению 1.98× относительно последовательной версии (0.18437 с).

 - Очень большие наборы данных (например, n > 100 000). При многократном увеличении n доля вычислений на одну точку растёт, а коммуникационные затраты (пересылка нескольких чисел на каждом шаге) становятся относительно малыми. При этом MPI позволяет распределить данные между узлами и ускорить обработку за счёт параллельной работы нескольких машин.

Когда гибридная схема не оправдана:

 - Использование двух и более MPI-процессов на задаче размера 2000 точек. Во всех конфигурациях с proc ≥ 2 время выполнения оказывается больше, чем у последовательной версии (ускорение < 1). Причина – высокие накладные расходы на коллективные операции (MPI_Gather, MPI_Bcast), вызываемые на каждом из ~2000 шагов алгоритма Джарвиса. Даже добавление OpenMP-потоков (конфигурация 2×4) не улучшает ситуацию – узким местом остаются коммуникации.

 - Любая задача, где размер данных мал (n < 10 000), а число шагов оболочки велико (h ≈ n). В таких условиях MPI только ухудшает производительность.