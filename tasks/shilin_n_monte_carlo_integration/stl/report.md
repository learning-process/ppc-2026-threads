# Многомерное интегрирование Монте-Карло — std::thread

- **Student:** Шилин Никита Дмитриевич, группа 3823Б1ПР1
- **Technology:** STL (`std::thread`)
- **Variant:** 12

---

## 1. Контекст

В STL-версии параллелизм управляется **вручную**: программа сама создаёт
пул потоков, распределяет работу между ними, объединяет частичные суммы
и завершает потоки через `join`. Это самая «инженерная» из четырёх
потоковых реализаций: всё, что в OMP/TBB прячется за директивами и
runtime-ом, здесь видно явно — диапазоны, локальные буферы, точка
синхронизации.

## 2. Постановка задачи

Та же, что в [`seq/report.md`](../seq/report.md) и
[`omp/report.md`](../omp/report.md). Тип задачи — `kSTL`, класс —
`ShilinNMonteCarloIntegrationSTL`.

## 3. Базовый алгоритм

Тот же замкнутый расчёт точки, что в OMP/TBB (без межитерационной
зависимости). Параллельная декомпозиция — **циклическая** по индексу:
поток `tid` обрабатывает индексы \(i = \mathrm{tid},\,\mathrm{tid}+T,\,\mathrm{tid}+2T,\dots\),
где \(T\) — общее число потоков. Это эквивалентно «странной» нарезке
диапазона `[0, n)`, но более устойчивой к разнородности итераций по
кэшу и хорошо ложится на квази-случайную последовательность Кронекера.

## 4. Схема распараллеливания

```49:93:tasks/shilin_n_monte_carlo_integration/stl/src/ops_stl.cpp
bool ShilinNMonteCarloIntegrationSTL::RunImpl() {
  auto dimensions = static_cast<int>(lower_bounds_.size());
  const std::vector<double> alpha = { /* sqrt(2)..sqrt(29) */ };

  // Honor PPC_NUM_THREADS so reports can build a real T-scaling table.
  auto num_threads = static_cast<unsigned int>(std::max(1, ppc::util::GetNumThreads()));

  std::vector<double> partial_sums(num_threads, 0.0);
  std::vector<std::thread> threads(num_threads);

  auto worker = [&](unsigned int tid) {
    std::vector<double> point(dimensions);
    for (int i = static_cast<int>(tid); i < num_points_; i += static_cast<int>(num_threads)) {
      for (int di = 0; di < dimensions; ++di) {
        double val = 0.5 + (static_cast<double>(i + 1) * alpha[di]);
        double current = val - std::floor(val);
        point[di] = lower_bounds_[di] + ((upper_bounds_[di] - lower_bounds_[di]) * current);
      }
      partial_sums[tid] += IntegrandFunction::Evaluate(func_type_, point);
    }
  };

  for (unsigned int ti = 0; ti < num_threads; ++ti) {
    threads[ti] = std::thread(worker, ti);
  }
  for (auto &th : threads) {
    th.join();
  }

  double sum = 0.0;
  for (unsigned int ti = 0; ti < num_threads; ++ti) {
    sum += partial_sums[ti];
  }
  ...
```

Ключевые места.

- **Источник числа потоков — `ppc::util::GetNumThreads()`**, то есть
  значение из переменной окружения `PPC_NUM_THREADS` (если она не
  задана — 1). Это согласовано с курсовым каркасом и позволяет
  построить осмысленную таблицу масштабирования по `T` (см. §8).
  Раньше использовался `std::thread::hardware_concurrency()` —
  переход на `GetNumThreads()` сделан в отдельном `[FIXED]` PR; до
  этой правки `PPC_NUM_THREADS` на STL-версию не влиял.
- **Локальные аккумуляторы.** Массив `partial_sums[num_threads]` —
  **поток-локальные** счётчики; каждый поток пишет только в свою ячейку
  `partial_sums[tid]`. Гонок по записи нет, потому что разные потоки
  обращаются к разным элементам.
- **False sharing — реальное узкое место.** `partial_sums[tid]` — это
  `double` (8 байт) подряд в `std::vector<double>`. На M4 Max длина
  кэш-линии 64 байта, то есть 8 соседних `double` лежат в одной линии.
  Когда 4–16 потоков в горячем цикле постоянно делают `+=` в свои
  `partial_sums[tid]`, кэш-линия с этими счётчиками постоянно мигрирует
  между ядрами через L2/L3-coherency-протокол. Это видно в таблицах §8:
  при `T = 4` эффективность падает до 59% (для сравнения OMP — 86%, TBB
  — 97%). «Правильное» исправление — `alignas(64)` обёртка на каждый
  счётчик, но в учебной версии оставлен простой массив; направление
  дальнейшей работы зафиксировано в выводах §9.
- **Создание потоков.** Цикл `threads[ti] = std::thread(worker, ti)`
  создаёт потоки **до** какого-либо `join`. Это сознательно сделано
  именно так, потому что распространённая ошибка (см. предупреждение
  в методичке курса о `example_threads/stl`) — вызывать `join`
  внутри цикла создания, что фактически сериализует работу.
- **`join`.** Отдельный цикл по `threads`; гарантирует, что главный
  поток дождётся всех рабочих перед чтением `partial_sums`. После
  объединения частичных сумм запись в `GetOutput()` идёт уже на одной
  нити — гонок нет.
- **Синхронизации между рабочими потоками внутри `worker`** нет: каждый
  поток работает в своём подмножестве индексов, в свой буфер `point` и в
  свою ячейку `partial_sums[tid]`. Это и есть «правильная» STL-версия
  Monte Carlo: ручная декомпозиция плюс детерминированная редукция в
  конце.

## 5. Детали реализации

- Файлы: [`stl/include/ops_stl.hpp`](include/ops_stl.hpp),
  [`stl/src/ops_stl.cpp`](src/ops_stl.cpp).
- Класс `ShilinNMonteCarloIntegrationSTL` возвращает `TypeOfTask::kSTL`.
- Число потоков читается через `ppc::util::GetNumThreads()`
  (`util/include/util.hpp`); если переменная не задана, fallback —
  1 поток. Это согласовано с курсовым каркасом и со стилем других
  STL-задач курса.

## 6. Проверка корректности

8 функциональных кейсов из общего набора курса
([`tests/functional/main.cpp`](../tests/functional/main.cpp)) — все
проходят. Эталон — аналитический интеграл; допуск ошибки тот же, что в
SEQ/OMP. Результат `sum` детерминирован при фиксированных
`num_threads` и `n` (потому что разбиение по индексу циклическое и
чтение детерминированное), но порядок суммирования зависит от
`num_threads` — поэтому тестовый порог
\(\max(10V/\sqrt n,\,10^{-2})\) выдержан с большим запасом.

## 7. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | Apple M4 Max, 16 ядер (12P + 4E) |
| RAM | 64 GiB |
| OS | macOS 26.3.1 |
| Компилятор | Apple clang 17.0.0 |
| Сборка | Release, `build-local` |
| Размер | \(n = 10^7\), \([0,1]^3\), `kSumSquares` |
| Источник `num_threads` | `ppc::util::GetNumThreads()` ↔ `PPC_NUM_THREADS` |
| Повторов | 5 запусков, медиана |

Команда замера для конкретного `T`:

```bash
export PPC_NUM_THREADS=$T OMP_NUM_THREADS=$T
./build-local/bin/ppc_perf_tests --gtest_filter='*shilin_n_monte_carlo_integration_stl*'
```

## 8. Результаты

Базовое `T_seq(task_run) = 0,030436 c`,
`T_seq(pipeline) = 0,030642 c` (см. [`seq/report.md`](../seq/report.md)).

### 8.1. Режим `task_run`

| `PPC_NUM_THREADS` | время, с | \(S\) | \(E\) |
| ---: | ---: | ---: | ---: |
| 1 | 0,030621 | 0,99 | 99,4% |
| 2 | 0,025104 | 1,21 | 60,6% |
| 4 | 0,012939 | 2,35 | 58,8% |
| 8 | 0,007340 | 4,15 | 51,8% |
| 16 | **0,006121** | **4,97** | 31,1% |

### 8.2. Режим `pipeline`

| `PPC_NUM_THREADS` | время, с | \(S\) | \(E\) |
| ---: | ---: | ---: | ---: |
| 1 | 0,030647 | 1,00 | 100,0% |
| 2 | 0,019838 | 1,54 | 77,2% |
| 4 | 0,015251 | 2,01 | 50,2% |
| 8 | 0,007081 | 4,33 | 54,1% |
| 16 | **0,006413** | **4,78** | 29,9% |

**Наблюдения.**

1. **Видимый false sharing.** При `T = 2..4` эффективность опускается до
   58–60%, тогда как у OMP/TBB при тех же `T` остаётся 86–97%. Причина —
   соседние `partial_sums[tid]` лежат в одной 64-байтовой кэш-линии, и
   её `MOESI`-перемещения между ядрами съедают часть параллелизма
   (см. §5). Это ожидаемая «цена» простой ручной редукции на массив без
   `alignas(64)`.
2. **При `T = 16` `S ≈ 5`** — STL заметно проигрывает OMP (`S=8.5`) и
   TBB (`S=11.6`) при том же физическом числе ядер. К false sharing
   добавляются: (а) накладные расходы создания и `join` 16 потоков на
   каждый прогон perf-теста, (б) отсутствие пула — потоки создаются
   заново при каждом `RunImpl`, (в) отсутствие закрепления потоков за
   ядрами — ОС свободно мигрирует их между P и E.
3. `task_run` и `pipeline` дают близкие, но не одинаковые цифры — это
   характерно для STL из-за переменных накладных расходов запуска
   потоков; разница укладывается в шум.

## 9. Выводы

- STL-версия **корректна** (8/8 функциональных тестов) и
  масштабируется по `PPC_NUM_THREADS` (`S = 5,0` при `T=16`).
- Однако по эффективности она **уступает** OMP и TBB: false sharing на
  `partial_sums[]` + накладные расходы создания 16 системных потоков на
  каждый `RunImpl` дают `E ≈ 31%` при `T=16`. Это типичная «цена»
  простого ручного параллелизма без оптимизаций.
- **Основное направление дальнейшей работы** — перейти на
  `alignas(64)`-обёртки на счётчики (устранит false sharing) и
  переиспользовать пул потоков (через `std::jthread` + `std::stop_token`
  или однажды созданные long-living `std::thread` с queue) — это можно
  поднять `E` до уровня OMP. В данной версии оставлен «учебный»
  вариант, чтобы было видно цену прямолинейного решения.
- Несмотря на это, STL-версия годится в качестве **инженерного эталона**
  ручного параллелизма и проходит все формальные требования курса.

## 10. Источники

- Документация курса PPC, репозиторий `ppc-2026-threads`.
- C++ reference: [`std::thread`](https://en.cppreference.com/w/cpp/thread/thread),
  [`std::thread::hardware_concurrency`](https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency).
- Microsoft Docs: [Thread safety of standard library types](https://learn.microsoft.com/en-us/cpp/standard-library/thread-safety-in-the-cpp-standard-library).

## 11. Чек-лист

- [x] Описано циклическое разбиение по индексу, локальные буферы,
      `partial_sums`.
- [x] Подчёркнуто, что `join` происходит **после** всех создания всех
      потоков (а не внутри цикла создания).
- [x] Зафиксирован переход с `hardware_concurrency()` на
      `ppc::util::GetNumThreads()` (`[FIXED]` PR), с описанием эффекта
      (теперь STL уважает `PPC_NUM_THREADS`).
- [x] Корректность сверена с SEQ при `T ∈ {1, 8, 16}`.
- [x] Таблицы по \(T \in \{1, 2, 4, 8, 16\}\) для `task_run` и
      `pipeline` (медиана 5).
- [x] Объяснён видимый false sharing на `partial_sums` и направление
      дальнейшей оптимизации (`alignas(64)`).
- [x] Этот файл — `stl/report.md`, лежит ровно в `tasks/<task>/stl/`.
