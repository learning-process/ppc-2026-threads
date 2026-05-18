# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3. — STL

- Student: Исхаков Дамир Айратович
- Technology: STL
- Variant: 25


## 1. Контекст

Версия STL реализует параллелизм с помощью стандартных потоков `std::thread`.
Изображение вручную разбивается на вертикальные полосы (blocks), каждый поток обрабатывает свой диапазон 
столбцов независимо. Это позволяет оценить производительность ручного управления потоками в сравнении с 
высокоуровневыми технологиями (OMP, TBB)


## 2. Постановка задачи

Постановка задачи ровно такая же как и в `seq/report.md` (пункт 2). включая **входные данные**, **выходные данные**, 
**ограничения** и **крайние случаи**


## 3. Базовый алгоритм

Последовательный алгоритм тот же, что и в SEQ: свёртка 3×3 с ядром Гаусса и зеркальным отражением границ
Асимптотика по времени O(N), по памяти O(N). Подробнее см. `seq/report.md` (пункт 3)


## 4. Схема распараллеливания

- **Декомпозиция данных:** - изображение делится на `actual_threads` вертикальных полос. Ширина каждой полосы 
вычисляется как `cols_per_thread = width / actual_threads`, остаток `remainder` распределяется по первым потокам.
- **Потоки:** - для каждой полосы создаётся отдельный `std::thread` с лямбда-функцией, обрабатывающей свои столбцы 
от `start_col` до `end_col`.
- **Локальные результаты** - каждый поток читает общую матрицу `matrix` и пишет в уникальные ячейки `result`, 
поэтому **синхронизация не требуется**
- **Join:** - все потоки создаются в цикле и сохраняются в `std::vector<std::thread>`, а затем последовательно 
вызывается `join()` для каждого. Это корректный паттерн: потоки работают параллельно, а главный поток ждёт их 
завершения.


## 5. Детали реализации

Файлы: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`

В сравнении с SEQ, OMP и TBB версия STL практически идентична:

* **ValidationImpl** и **PreProcessingImpl**/**PostProcessingImpl** – без изменений относительно SEQ, OMP и TBB.
* **RunImpl** изменения относительно SEQ, OMP и TBB: 

1. Определяется число потоков `num_threads = ppc::util::GetNumThreads()` и ограничивается шириной изображения: 
`actual_threads = std::min(num_threads, width)`.
2. Вычисляются границы полос для каждого потока и в цикле создаются `std::thread`:

```cpp
for (int thread_id = 0; thread_id < actual_threads; ++thread_id) {
    int end_col = start_col + cols_per_thread + (thread_id < remainder ? 1 : 0);
    threads.emplace_back([&, start_col, end_col]() {
      for (int horizontal_band = start_col; horizontal_band < end_col; ++horizontal_band) {
        for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
          int sum = 0;

          sum += kGaussKernel[0][0] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band - 1, vertical_band - 1, width, height);
          sum += kGaussKernel[0][1] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band, vertical_band - 1, width, height);
          sum += kGaussKernel[0][2] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band + 1, vertical_band - 1, width, height);

          sum += kGaussKernel[1][0] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band - 1, vertical_band, width, height);
          sum += kGaussKernel[1][1] * 
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band, vertical_band, width, height);
          sum += kGaussKernel[1][2] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band + 1, vertical_band, width, height);

          sum += kGaussKernel[2][0] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band - 1, vertical_band + 1, width, height);
          sum += kGaussKernel[2][1] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band, vertical_band + 1, width, height);
          sum += kGaussKernel[2][2] *
                 IskhakovDGetPixelMirrorSTL(matrix, horizontal_band + 1, vertical_band + 1, width, height);

          result[(vertical_band * width) + horizontal_band] = static_cast<uint8_t>(sum / kDivConst);
        }
      }
    });
    start_col = end_col;
  }
```

Функция **IskhakovDGetPixelMirrorStl** осталась без изменений, поменялось только название реализации в конце 
(вместо Seq, Omp или Tbb)


## 6. Проверка корректности

Набор тестов для STL идентичен набору тестов для SEQ, и также проходят все функциональные тесты 
(см пункт 5 `seq/report.md`)


## 7. Экспериментальная среда

- **CPU:** 12th Gen Intel(R) Core(TM) i5-12500H (12 физических ядер, 16 логических потоков, 
базовая частота 400 МГц – 4.5 ГГц)
- **RAM:** 16 ГБ
- **OS:** Linux Mint 22.3 (Zena)
- **Компилятор:** GCC 15.2.0
- **Сборка:** CMake 3.28.3, Release
- **Команда запуска функциональных тестов:** 
`./build/bin/ppc_func_tests --gtest_filter="*IskhakovDVerticalGaussFilterFuncTests*stl*"`
- **Команда запуска тестов производительности:** 
  ```bash
       PPC_NUM_THREADS=<число потоков> ./build/bin/ppc_perf_tests 
       --gtest_filter="*IskhakovDVerticalGaussFilterPerfTests*stl*"
  ```

## 8. Результаты

Измерение производительности проводилось на изображении размером 8192×8192 пикселей (~67 млн пикселей).  
Базовое время последовательной версии (SEQ) – 1.79 с (ускорение = 1, эффективность = 100%).

| Число потоков | Время выполнения, с | Ускорение (Speedup) | Эффективность (Efficiency) |
| ------------- | ------------------- | ------------------- | -------------------------- |
| 1 (SEQ) | 1.79 | 1.00 | 100.0% |
| 2 | 1.05 | 1.70 | 85.0% |
| 4 | 0.74 | 2.42 | 60.5% |
| 8 | 0.79 | 2.27 | 28.4% |
| 12 | 0.76 | 2.36 | 19.7% |
| 16 | 0.80 | 2.24 | 14.0% |

Коэффициент ускорения: S = T_seq / T_p. Эффективность: E = S / p × 100%.

Самым времязатратным фрагментом остаётся двойной цикл for.
Версия на std::thread демонстрирует уверенный рост ускорения до 4 потоков (2.42×), после чего время практически 
не меняется с увеличением числа потоков
Это объясняется накладными расходами на создание и join большего числа потоков
По итогу реализация STL показала низкую максимальную производительность (2.42× против 2.67× у OMP и 2.91× у TBB) 
и заметное падение эффективности при большом числе потоков


## 9. Выводы

- Ручная реализация на std::thread обеспечила максимальное ускорение **2.42×** (на 4 потоках), 
немного уступив OMP (2.67×) и TBB (2.91×).
- Наилучший баланс достигается при 4 потоках: ускорение 2.42, эффективность 60.5%.
- При увеличении числа потоков свыше 4 производительность практически не изменяется, а на 16 потоках даже снижается 
из-за накладных расходов на создание/join потоков
- std::thread предоставляет полный контроль над разбиением данных и жизненным циклом потоков, но требует ручного 
кодирования и более чувствителен к дисбалансу нагрузки.
- В условиях данной задачи высокоуровневые технологии (OMP, TBB) оказались более эффективными, особенно на большом 
количестве потоков, благодаря встроенным механизмам балансировки и управления пулом потоков