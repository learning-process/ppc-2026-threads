# Умножение плотных матриц. Элементы типа `double`. Алгоритм Штрассена — TBB

- Student: Мухаммадхон Исрам
- Technology: TBB
- Variant: 3

## 1. Контекст

У Штрассена хорошо видны независимые ветви на верхнем уровне,
поэтому oneTBB здесь выглядит уместно. `parallel_invoke` позволяет
запустить все семь произведений параллельно без ручного управления потоками.

## 2. Постановка задачи

Та же постановка, что и в `SEQ`: `C = A * B`, допуск `1e-9`.

## 3. Базовый алгоритм

Математическая часть не менялась. `kCutoff = 64`, `kBlockSize = 64`.
Внутренняя часть — последовательный `StrassenSeq` через `std::function impl`.

## 4. Схема распараллеливания

Декомпозиция: верхний уровень `M1..M7` — семь независимых задач.

Примитив: `oneapi::tbb::parallel_invoke`. Каждая ветвь оформлена
как отдельная лямбда. `parallel_for`, `blocked_range` и `partitioner`
не использовались — задачи дискретные, а не диапазон.

Контроль потоков: число рабочих потоков ограничивается через
`global_control::max_allowed_parallelism` по `ppc::util::GetNumThreads()`.

Синхронизация: `parallel_invoke` сам дожидается всех задач;
каждая ветвь пишет в свой буфер `m1..m7`, гонок нет.

## 5. Детали реализации

Файлы:

- `tbb/include/ops_tbb.hpp`
- `tbb/src/ops_tbb.cpp`

По сравнению с `SEQ`:

- добавлен `StrassenTopTbb(...)` с `parallel_invoke`;
- внутренняя часть оставлена в `StrassenSeq(...)`;
- рекурсия через `std::function impl` — без NOLINT.

## 6. Проверка корректности

```bash
cd build/bin
mpiexec -n 1 ./ppc_func_tests --gtest_filter=*muhammadkhon_i_stressen_alg_tbb*
```

Все функциональные тесты пройдены локально.

## 7. Экспериментальная среда

- CPU: Intel Core i5-10300H, 4 ядра / 8 потоков;
- RAM: 8 GB;
- OS: Fedora Linux (виртуальная машина);
- compiler: GCC, build type: `Release`.

```bash
mpiexec -n 1 ./ppc_perf_tests --gtest_filter=*muhammadkhon_i_stressen_alg_tbb*
```

## 8. Результаты

| Mode | Workers | Metric | Time, s | Speedup | Efficiency |
| --- | ---: | --- | ---: | ---: | ---: |
| seq | 1 | task_run | 1.124371 | 1.00 | N/A |
| tbb | 8 | task_run | 0.398617 | 2.82 | 35.2% |

## 9. Выводы

`TBB` показал лучшее время среди потоковых backend-ов. `parallel_invoke`
хорошо подходит для задачи с ограниченным числом независимых крупных ветвей.
Эффективность 35.2% — ожидаемо: семь задач на восемь потоков,
внутренняя часть последовательная.