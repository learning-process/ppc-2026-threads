# Повышение контраста полутонового изображения посредством линейной растяжки гистограммы — TBB

- Student: Гутянский Алексей Сергеевич, 3823Б1ФИ3
- Technology: TBB
- Variant: 28

## 1. Контекст

Библиотека Threading Building Blocks (TBB) реализует параллелизм на основе задач, что позволяет динамически
балансировать нагрузку.Для рассматриваемой задачи (массив байт, независимая обработка каждого элемента) подойдет распараллеливание
циклов через `parallel_for` и `parallel_reduce`.

## 2. Постановка задачи

(См. `seq/report.md`, раздел 2.)

## 3. Базовый алгоритм

Алгоритм идентичен SEQ, но для нахождения минимума/максимума используется `parallel_reduce`, а для преобразования – `parallel_for`.

## 4. Схема распараллеливания

**Поиск границ:**

- `tbb::parallel_reduce` над диапазоном `[0, sz)`.
- Тело функции – лямбда, обновляющая локальные `min`/`max`.
- Функция слияния – `std::min` для первого компонента, `std::max` для второго.

**Преобразование:**

- `tbb::parallel_for` с `tbb::blocked_range<size_t>(0, sz)`.
- В лямбде обрабатывается непрерывный поддиапазон индексов.
- Для случая `delta == 0` – аналогичный `parallel_for` с копированием.

Разделение работы между потоками в TBB динамическое. Для разбиения области итерирования установлен grain_size,
равный количеству итераций деленному на количество потоков, что позволяет приблизиться к статическому планированию из OpenMP.

## 5. Детали реализации

Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`

```cpp
// Поиск минимального и максимального элементов
auto [lower_bound, upper_bound] = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, sz, chunk_sz), std::make_pair(static_cast<uint8_t>(255), static_cast<uint8_t>(0)),
      [&input](const auto &range, auto init) {
    auto [local_min, local_max] = init;
    for (size_t i = range.begin(); i != range.end(); ++i) {
      local_min = std::min(local_min, input[i]);
      local_max = std::max(local_max, input[i]);
    }
    return std::make_pair(local_min, local_max);
  }, [](auto a, auto b) { return std::make_pair(std::min(a.first, b.first), std::max(a.second, b.second)); });

// Преобразование элементов
tbb::parallel_for(tbb::blocked_range<size_t>(0, sz, chunk_sz), [&](const auto &range) {
      for (auto idx = range.begin(); idx != range.end(); ++idx) {
        uint16_t old_value = input[idx];
        uint16_t new_value = (kMaxUint8 * (old_value - lower_bound)) / delta;
        output[idx] = static_cast<uint8_t>(new_value);
      }
    });

```

## 6. Проверка корректности

Функциональные тесты проходят для всех наборов данных. Результаты совпадают с SEQ.

## 7. Экспериментальная среда

Информацию о среде можно увидеть в (`seq/report.md`,раздел 6.)

TBB версия: 2022.0.0

Команда запуска:

```powershell
$env:PPC_NUM_THREADS=<кол-во потоков>
$env:OMP_NUM_THREADS=<кол-во потоков>

./build/bin/ppc_perf_tests --gtest_filter=RunModeTests/GutyanskyARunPerfTestsImgContrastIncr.RunPerfModes/*
```

## 8. Результаты (N = 100 000 000)

| Потоков | Время, мс | Ускорение | Эффективность |
| ------- | --------- | --------- | ------------- |
| 1       | 1092      | 0.75      | 75%           |
| 2       | 606       | 1.36      | 68%           |
| 4       | 367       | 2.24      | 56%           |
| 6       | 297       | 2.77      | 46%           |
| 12      | 226       | 3.65      | 30%           |
| 16      | 220       | 3.75      | 23%           |

Ожидалось, что TBB с динамическим планированием даст результат, немного отстающий от того, который выдал OpenMP.
Но поддержка omp в компиляторе MSVC ограничена, потому результат TBB в среднем немного лучше чем у OpenMP.

## 9. Выводы

TBB удобен для задач, где важна масштабируемость на неоднородных данных. Для данной равномерной нагрузки разница
с OpenMP незначительна. Рекомендуется, когда требуется избежать привязки к компилятору
или использовать продвинутые планировщики задач.
