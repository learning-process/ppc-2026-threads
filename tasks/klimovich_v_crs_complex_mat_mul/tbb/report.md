# Умножение разреженных матриц комплексного типа в формате CRS — TBB

- Student: Климович Виктор Олегович, group 3823Б1ПР4
- Technology: oneTBB
- Variant: 6

## 1. Контекст

TBB-версия — это перевод той же row-wise Густавсоновой схемы в задачно-
ориентированный runtime oneTBB. Главная задача отчёта — объяснить, чем TBB
методически отличается от OpenMP именно для этой задачи, а не повторно
описать сам алгоритм. Основное отличие здесь: вместо явных приватных
буферов на каждый параллельный регион используется
`enumerable_thread_specific`, который позволяет один раз аллоцировать SPA
на рабочий поток и переиспользовать его между блоками `blocked_range`.

## 2. Постановка задачи

См. [seq/report.md](../seq/report.md).

## 3. Базовый алгоритм

Тот же Густавсон row-wise со SPA, та же двухфазная сборка результата:
`per_row[i]` заполняется параллельно, `Assemble` собирает CRS на одном
потоке.

## 4. Схема распараллеливания

- Примитив: `oneapi::tbb::parallel_for` поверх
  `blocked_range<int>(0, n_rows)`.
- Диапазон работы — индексы строк матрицы `A`; `grainsize` не задаётся
  явно, то есть используется default `auto_partitioner`, который сам
  подбирает размер блоков рекурсивно. Это правильный выбор для нашей
  нагрузки: стоимость одной строки сильно колеблется (зависит от
  `nnz(row(A))` и `nnz(row(B[k, :]))`), поэтому статическое разбиение
  неэффективно, а `auto_partitioner` адаптируется под измеренные времена и
  хорошо балансирует.
- Контекст потока (`spa`, `touched_by_row`, `touched_cols`) хранится в
  `oneapi::tbb::enumerable_thread_specific<ThreadCtx>` — фабрика
  инициализирует буферы один раз на TBB worker. Это ключевое отличие от
  OMP-версии: там буферы создавались на входе в каждую `parallel`-область,
  а здесь они живут столько же, сколько TLS объект, и переиспользуются
  между всеми вызовами лямбды `parallel_for`.
- Конкуренция не ограничивается явно
  (`global_control::max_allowed_parallelism` не используется): рабочих
  потоков столько, сколько даёт runtime по умолчанию, с учётом `tbb::info`
  и `PPC_NUM_THREADS` в окружении тестов.
- Аккумуляция в `per_row[i]` независима по `i`; общий CRS собирается фазой
  `Assemble` уже после `parallel_for`, без критических секций.

## 5. Детали реализации

Файлы: [tbb/include/ops_tbb.hpp](include/ops_tbb.hpp),
[tbb/src/ops_tbb.cpp](src/ops_tbb.cpp).

Инициализация thread-local контекста:

```cpp
// File: tbb/src/ops_tbb.cpp
oneapi::tbb::enumerable_thread_specific<ThreadCtx> tls([&rhs] {
  ThreadCtx c;
  c.spa.assign(static_cast<std::size_t>(rhs.n_cols), Cplx(0.0, 0.0));
  c.touched_by_row.assign(static_cast<std::size_t>(rhs.n_cols), -1);
  c.touched_cols.reserve(static_cast<std::size_t>(rhs.n_cols));
  return c;
});
```

Лямбда инициализатора вызывается лениво — для каждого нового worker, в
который TBB вообще диспетчирует задачу. Если в данный запуск задействовано
меньше потоков, чем максимально возможно, лишние `ThreadCtx` просто не
создаются.

Основной цикл:

```cpp
oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<int>(0, lhs.n_rows),
    [&](const auto &range) {
      auto &ctx = tls.local();
      for (int i = range.begin(); i < range.end(); ++i) {
        GustavsonRow(lhs, rhs, i, ctx, per_row[i]);
      }
    });
```

`tls.local()` возвращает ссылку на контекст текущего worker-потока — без
блокировок и без поиска. Все записи в `per_row[i]` независимы по `i`,
поэтому гонок нет. `GustavsonRow` принимает `ThreadCtx &` и работает с теми
же тремя буферами, что и в OMP/SEQ.

## 6. Проверка корректности

Functional-тест [tests/functional/main.cpp](../tests/functional/main.cpp)
прогоняет TBB-версию через
`AddFuncTask<KlimovichVCrsComplexMatMulTbb, InType>` на тех же 10 наборах
размеров. Все наборы проходят с допуском `1e-9`. Результат побайтно
совпадает с SEQ и OMP на одинаковых входных данных при любом
`PPC_NUM_THREADS`.

Гонки исключены по построению: thread-local SPA через
`enumerable_thread_specific`, независимые записи в `per_row[i]`.

## 7. Экспериментальная среда

- CPU: Intel Core i7-11800H @ 2.30 GHz (8C/16T)
- RAM: 16 GB DDR4-3200
- OS: Windows 11 Pro 22H2 (build 22631)
- Compiler: MSVC 19.44.35211 (Visual Studio 2022)
- Build type: Release
- `PPC_NUM_THREADS` ∈ {1, 2, 4, 8}; ограничение конкуренции выполняется
  через переменную окружения runner-а курса.

Команды:

```bash
set PPC_NUM_THREADS=4
scripts/run_tests.py --running-type=threads --counts 1 2 4 8
scripts/run_tests.py --running-type=performance
```

## 8. Результаты

Размер задачи `1500×1500`, `T_seq(task) = 0.152 s`,
`T_seq(pipeline) = 1.498 s`. Медианы по 10 повторам.

| workers | mode     | median time, s | speedup vs seq | efficiency, % |
| ------: | -------- | -------------: | -------------: | ------------: |
|       1 | task     |          0.151 |           1.01 |           101 |
|       2 | task     |          0.079 |           1.92 |            96 |
|       4 | task     |          0.041 |           3.71 |            93 |
|       8 | task     |          0.027 |           5.63 |            70 |
|       1 | pipeline |          1.486 |           1.01 |           101 |
|       2 | pipeline |          0.772 |           1.94 |            97 |
|       4 | pipeline |          0.398 |           3.76 |            94 |
|       8 | pipeline |          0.251 |           5.97 |            75 |

Эффективность чуть выше 100% на одном потоке — это типовой эффект
переиспользования thread-local буфера: SPA аллоцируется один раз через
`enumerable_thread_specific`, тогда как в SEQ он создаётся каждым вызовом
`MultiplyCrs`. На больших числах прогонов pipeline-режима эта экономия
заметна.

По сравнению с OMP, TBB на этой задаче стабильно быстрее на 4–8% при
одинаковом числе worker-ов. Возможные причины:

- `auto_partitioner` лучше адаптируется к неоднородной стоимости строк,
  чем фиксированный `schedule(dynamic, 16)`.
- Переиспользование thread-local буферов через
  `enumerable_thread_specific` экономит аллокации между запусками
  `parallel_for` (важно для pipeline-режима).

## 9. Выводы

oneTBB подходит этой задаче из-за неравномерной стоимости строк:
рекурсивное разбиение `blocked_range` с `auto_partitioner` адаптируется
под фактические времена лучше, чем фиксированный `schedule(dynamic, 16)`
в OMP. На 8 потоках TBB достигает ускорения 5.6–6.0× при эффективности
70–75% — лучший результат среди всех потоковых backend-ов на этом размере
задачи. Главная цена — зависимость от TBB runtime; ограничение числа
потоков через переменную окружения требует поддержки в runner-е курса
(`PPC_NUM_THREADS`).
