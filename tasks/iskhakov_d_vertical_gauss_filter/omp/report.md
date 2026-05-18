# Линейная фильтрация изображений(вертикальное разбиение).Ядро Гаусса 3x3. — OMP

- Student: Исхаков Дамир Айратович
- Technology: OMP
- Variant: 25

## 1. Контекст

Из SEQ в OMP переносится вся реализация, только вместо "формального" разбиения на полосы реализации SEQ в OMP
реализуется полноценное разбиение на полосы относительно числа потоков, благодаря чему задача выполняется параллельно

## 2. Постановка задачи

Постановка задачи ровно такая же как и в `seq/report.md` (пункт 2). включая **входные данные**, **выходные данные**,
**ограничения** и **крайние случаи**

## 3. Базовый алгоритм

Последовательный алгоритм тот же, что и в SEQ: свёртка 3×3 с ядром Гаусса и зеркальным отражением границ
Асимптотика по времени O(N), по памяти O(N). Подробнее см. `seq/report.md` (пункт 3)

## 4. Схема распараллеливания

**Область распараллеливания** - в данной реализации распараллеливается сугубо внешний цикл по столбцам

```cpp
#pragma omp parallel for default(none) shared(matrix, result, width, height, kGaussKernel, kDivConst)
  for (int horizontal_band = 0; horizontal_band < width; ++horizontal_band)
```

**Переменные:**

1. `shared: matrix, result, width, height, kGaussKernel, kDivConst` - общие для всех потоков
2. `private: horizontal_band, vertical_band, sum` - локальные в каждом потоке

**Синхронизация (Где нужен reduction / atomic / critical)** - в данной реализации синхронизация не требуется,
т.к каждый поток работает со своим набором данных и изменяет сугубо свои ячейки в векторе результата (изменение
и определение своих ячеек происходит за счет обращения к их индексам horizontal_band и vertical_band)

**Планирование (schedule)** - в реализации используется планировщик по умолчанию static, благодаря чему потоки сразу
получают их полосы и необходимые данные для выполнения задачи, без дополнительных изменений во время работы

**Число потоков** - задается при запуске через параметр `PPC_NUM_THREADS` (или `OMP_NUM_THREADS`) и передаётся через
`omp_set_num_threads`

## 5. Детали реализации

Файлы: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`

В сравнении с SEQ версия OMP практически идентична:

- **ValidationImpl** и **PreProcessingImpl**/**PostProcessingImpl** – без изменений относительно SEQ.
- **RunImpl** изменения относительно SEQ:

1. Определение числа потоков через `omp_set_num_threads(num_threads);`
2. Была добавлена область распараллеливания:

```cpp
#pragma omp parallel for default(none) shared(matrix, result, width, height, kGaussKernel, kDivConst)
  for (int horizontal_band = 0; horizontal_band < width; ++horizontal_band) {
    for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
      int sum = 0;

      sum += kGaussKernel[0][0] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band - 1, vertical_band - 1, width, height);
      sum += kGaussKernel[0][1] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band, vertical_band - 1, width, height);
      sum += kGaussKernel[0][2] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band + 1, vertical_band - 1, width, height);

      sum += kGaussKernel[1][0] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band - 1, vertical_band, width, height);
      sum += kGaussKernel[1][1] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band, vertical_band, width, height);
      sum += kGaussKernel[1][2] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band + 1, vertical_band, width, height);

      sum += kGaussKernel[2][0] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band - 1, vertical_band + 1, width, height);
      sum += kGaussKernel[2][1] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band, vertical_band + 1, width, height);
      sum += kGaussKernel[2][2] *
             IskhakovDGetPixelMirrorOmp(matrix, horizontal_band + 1, vertical_band + 1, width, height);

      result[(vertical_band * width) + horizontal_band] = static_cast<uint8_t>(sum / kDivConst);
    }
  }
```

Функция **IskhakovDGetPixelMirrorOmp** осталась без изменений, только название поменялось с
**IskhakovDGetPixelMirrorSeq** на **IskhakovDGetPixelMirrorOmp**

## 6. Проверка корректности

Набор тестов для OMP идентичен набору тестов для SEQ, и также проходят все функциональные тесты
(см пункт 5 `seq/report.md`)

## 7. Экспериментальная среда

- **CPU:** 12th Gen Intel(R) Core(TM) i5-12500H (12 физических ядер, 16 логических потоков,
базовая частота 400 МГц – 4.5 ГГц)
- **RAM:** 16 ГБ
- **OS:** Linux Mint 22.3 (Zena)
- **Компилятор:** GCC 15.2.0
- **Сборка:** CMake 3.28.3, Release
- **Команда запуска функциональных тестов:**
`./build/bin/ppc_func_tests --gtest_filter="*IskhakovDVerticalGaussFilterFuncTests*omp*"`
- **Команда запуска тестов производительности:**

  ```bash
  PPC_NUM_THREADS=<число потоков> ./build/bin/ppc_perf_tests \
    --gtest_filter="*IskhakovDVerticalGaussFilterPerfTests*omp*"
  ```

## 8. Результаты

Измерение производительности проводилось на изображении размером 8192×8192 пикселей (~67 млн пикселей)
Базовое время последовательной версии (SEQ) – 1.79 с (ускорение = 1, эффективность = 100%).

| Число потоков | Время выполнения, с | Ускорение (Speedup) | Эффективность (Efficiency) |
| ------------- | ------------------- | ------------------- | -------------------------- |
| 1 (SEQ) | 1.79 | 1.00 | 100.0% |
| 2 | 1.02 | 1.75 | 87.7% |
| 4 | 0.67 | 2.67 | 66.8% |
| 8 | 0.73 | 2.45 | 30.6% |
| 12 | 0.69 | 2.59 | 21.6% |
| 16 | 0.68 | 2.63 | 16.5% |

Коэффициент ускорения рассчитывается как S = T_seq / T_p, коэффициент эффективности — как E = S / p × 100%.

Самым времязатратным фрагментом, как и в SEQ, остаётся двойной цикл for с обращением к пикселям.
Начиная с 8 потоков время выполнения хоть и снижается, но эффективность резко падает из-за роста накладных расходов
OpenMP на распределение итераций и синхронизацию. Это говорит о том, что использование более 4 потоков для данной
задачи малоэффективно: основное время начинает уходить на пересылку данных и организацию параллельной работы.

## 9. Выводы

Выигрышем в данной ситуации можно считать использование 4 потоков, так как это решение даёт высокую скорость выполнения
(ускорение 2.67) при сохранении приемлемой эффективности (66.8%) и с минимальными накладными расходами
