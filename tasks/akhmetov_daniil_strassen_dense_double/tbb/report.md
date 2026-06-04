# Умножение плотных матриц. Элементы типа `double`. Алгоритм Штрассена

- Студент: Ахметов Даниил Данисович, группа 3823Б1ПР2
- Технология: TBB
- Вариант: 3

## 1. Введение

У Штрассена на верхнем уровне есть семь независимых ветвей `M1..M7`,
поэтому oneTBB здесь выглядит естественно. В этой версии верхний уровень
распараллеливается через `parallel_invoke`, внутренняя часть остаётся
последовательной с блочным умножением.

Общий контекст - в [report.md](../report.md).

## 2. Постановка задачи

Та же постановка, что и в `SEQ`: `C = A * B`, плотные матрицы `double`,
допуск `1e-7`. Подробнее - в [seq/report.md](../seq/report.md).

## 3. Описание алгоритма (базового/последовательного)

Внутренняя часть - `StrassenSeqImpl` с блочным умножением
(`kCutoff = 256`, `kBlockSize = 64`). Математика Штрассена та же,
что в `SEQ`.

## 4. Схема распараллеливания

**Декомпозиция:** верхний уровень `M1..M7` - семь независимых задач.

**Примитив:** `oneapi::tbb::parallel_invoke` - каждая ветвь
оформлена отдельной лямбдой.

**Планирование:**

- `parallel_for` не используется - задачи дискретные, не диапазон;
- `blocked_range`, `grainsize`, `partitioner` не применяются.

**Контроль потоков:**

- `global_control::max_allowed_parallelism`
  по `ppc::util::GetNumThreads()`.

**Синхронизация:**

- `parallel_invoke` сам дожидается всех задач;
- каждая ветвь пишет в свой буфер `m1..m7`, гонок нет.

## 5. Детали реализации

**Файлы:**

- `tbb/include/ops_tbb.hpp`
- `tbb/src/ops_tbb.cpp`

**Ключевые функции:**

- `StrassenTopTbb(...)` - верхний уровень с `parallel_invoke`;
- `StrassenSeqImpl(...)` - внутренняя последовательная часть;
- `NaiveMulBlocked(...)` - блочное умножение на малых размерах.

Матрицы адресуются через указатели и `stride`.

## 6. Экспериментальная установка

- CPU: AMD Ryzen 7 2700, 8 ядер / 16 потоков;
- RAM: 16 GB;
- OS: Windows 10 x64;
- IDE: Visual Studio Code;
- compiler: MSVC (сборка через CMake);
- build type: `Release`;
- `PPC_NUM_THREADS=4`.

**Команды:**

```powershell
cd build\bin
$env:PPC_NUM_THREADS = "4"
.\ppc_func_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*_tbb_*
.\ppc_perf_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double_tbb_enabled*
```

## 7. Результаты и обсуждение

### 7.1 Корректность

- размеры `64`, `128`, `256`;
- сравнение с наивным эталоном;
- локально все тесты пройдены.

### 7.2 Производительность

| Режим | Число потоков | Время, с | Ускорение | Эффективность |
| --- | ---: | ---: | ---: | ---: |
| seq | 1 | 0.938515 | 1.00 | N/A |
| tbb | 4 | 0.103911 | 9.03 | 225.8% |

`TBB` показал лучшее локальное время среди параллельных backend-ов.
Эффективность выше 100% связана с другой архитектурой верхнего уровня,
а не только с добавлением потоков к схеме `SEQ`.

## 8. Заключение

`parallel_invoke` хорошо подошёл для верхнего уровня Штрассена.
`TBB` - один из самых быстрых thread-based вариантов на локальных замерах.

## 9. Источники

1. [oneTBB Documentation](https://uxlfoundation.github.io/oneTBB/)
2. [Parallel Programming Course](https://learning-process.github.io/parallel_programming_course/ru/common_information/threading_tasks.html)

## Приложение (опционально)

```cpp
oneapi::tbb::parallel_invoke(
    [&] { /* M1 */ }, [&] { /* M2 */ }, /* ... M7 ... */
    [&] { CombineQuadrants(...); });
```
