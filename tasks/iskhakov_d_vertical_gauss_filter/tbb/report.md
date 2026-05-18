# Линейная фильтрация изображений(вертикальное разбиение).Ядро Гаусса 3x3. — TBB

- Student: Исхаков Дамир Айратович
- Technology: TBB
- Variant: 25


## 1. Контекст

Версия TBB переносит SEQ(последовательный алгоритм) фильтрации на технологию Intel oneTBB. Параллелизм достигается 
с помощью примитива `tbb::parallel_for`, который автоматически разбивает внешний цикл по столбцам изображения на 
независимые задачи и распределяет их между рабочими потоками


## 2. Постановка задачи

Постановка задачи ровно такая же как и в `seq/report.md` (пункт 2). включая **входные данные**, **выходные данные**, 
**ограничения** и **крайние случаи**


## 3. Базовый алгоритм
Последовательный алгоритм тот же, что и в SEQ: свёртка 3×3 с ядром Гаусса и зеркальным отражением границ
Асимптотика по времени O(N), по памяти O(N). Подробнее см. `seq/report.md` (пункт 3)


## 4. Схема распараллеливания

- **Примитив** - Для данной реализации был выбран `parallel_for`, который автоматически получает весь размер ширины 
`[0, width)`, а также тело цикла, после чего идет параллельный подсчет значений пикселей для каждой полосы 
`horizontal_band`
- **Диапазон** - интервал `[0, width)` (индексы пикселей от 0 до максимальной ширины изображения) передаваемый 
в примитив `parallel_for`
- **Grainsize** явно не задан, используется автоматическое разбиение
- **Partitioner** - по умолчанию применяется `auto_partitioner`, который рекурсивно делит диапазон, пока блоки 
не станут достаточно мелкими, благодаря чему нет необходимости делать это вручную
- **Ограничение конкуренции** задается через `PPC_NUM_THREADS <число потоков>` 


## 5. Детали реализации

Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`

В сравнении с SEQ и OMP версия TBB практически идентична:

* **ValidationImpl** и **PreProcessingImpl**/**PostProcessingImpl** – без изменений относительно SEQ и OMP.
* **RunImpl** изменения относительно SEQ и OMP: 

1. Вместо двойного цикла for используется `tbb::parallel_for(0, width, [&](int horizontal_band)`:

```cpp
tbb::parallel_for(0, width, [&](int horizontal_band) {
    for (int vertical_band = 0; vertical_band < height; ++vertical_band) {
      int sum = 0;

      sum += kGaussKernel[0][0] *
             IskhakovDGetPixelMirrorTbb(matrix, horizontal_band - 1, vertical_band - 1, width, height);
      sum += kGaussKernel[0][1] * 
              IskhakovDGetPixelMirrorTbb(matrix, horizontal_band, vertical_band - 1, width, height);
      sum += kGaussKernel[0][2] *
             IskhakovDGetPixelMirrorTbb(matrix, horizontal_band + 1, vertical_band - 1, width, height);

      sum += kGaussKernel[1][0] * 
              IskhakovDGetPixelMirrorTbb(matrix, horizontal_band - 1, vertical_band, width, height);
      sum += kGaussKernel[1][1] * 
              IskhakovDGetPixelMirrorTbb(matrix, horizontal_band, vertical_band, width, height);
      sum += kGaussKernel[1][2] * 
              IskhakovDGetPixelMirrorTbb(matrix, horizontal_band + 1, vertical_band, width, height);

      sum += kGaussKernel[2][0] *
             IskhakovDGetPixelMirrorTbb(matrix, horizontal_band - 1, vertical_band + 1, width, height);
      sum += kGaussKernel[2][1] * 
              IskhakovDGetPixelMirrorTbb(matrix, horizontal_band, vertical_band + 1, width, height);
      sum += kGaussKernel[2][2] *
             IskhakovDGetPixelMirrorTbb(matrix, horizontal_band + 1, vertical_band + 1, width, height);

      result[(vertical_band * width) + horizontal_band] = static_cast<uint8_t>(sum / kDivConst);
    }
  });
```

Функция **IskhakovDGetPixelMirrorTbb** осталась без изменений, поменялось только название реализации в конце 
(вместо Seq или Omp)


## 6. Проверка корректности

Набор тестов для TBB идентичен набору тестов для SEQ, и также проходят все функциональные тесты 
(см пункт 5 `seq/report.md`)


## 7. Экспериментальная среда

- **CPU:** 12th Gen Intel(R) Core(TM) i5-12500H (12 физических ядер, 16 логических потоков, 
базовая частота 400 МГц – 4.5 ГГц)
- **RAM:** 16 ГБ
- **OS:** Linux Mint 22.3 (Zena)
- **Компилятор:** GCC 15.2.0
- **Сборка:** CMake 3.28.3, Release
- **Команда запуска функциональных тестов:** 
`./build/bin/ppc_func_tests --gtest_filter="*IskhakovDVerticalGaussFilterFuncTests*tbb*"`
- **Команда запуска тестов производительности:** 

  ```bash
       PPC_NUM_THREADS=<число потоков> ./build/bin/ppc_perf_tests 
       --gtest_filter="*IskhakovDVerticalGaussFilterPerfTests*tbb*"
  ```



## 8. Результаты

Измерение производительности проводилось на изображении размером 8192×8192 пикселей (~67 млн пикселей)
Базовое время последовательной версии (SEQ) – 1.79 с (ускорение = 1, эффективность = 100%).

| Число потоков | Время выполнения, с | Ускорение (Speedup) | Эффективность (Efficiency) |
| ------------- | ------------------- | ------------------- | -------------------------- |
| 1 (SEQ) | 1.79 | 1.00 | 100.0% |
| 2 | 1.07 | 1.68 | 83.9% |
| 4 | 0.65 | 2.74 | 68.5% |
| 8 | 0.71 | 2.51 | 31.4% |
| 12 | 0.63 | 2.84 | 23.7% |
| 16 | 0.62 | 2.91 | 18.2% |

Коэффициент ускорения рассчитывается как S = T_seq / T_p, коэффициент эффективности — как E = S / p × 100%.

Самым времязатратным фрагментом остаётся двойной цикл for с обращением к пикселям.  
TBB показал небольшое преимущество в максимальном ускорении (2.91×) по сравнению с OMP (2.67×), при этом лучше 
сгладил провал на 8 потоках (2.51× против 2.45× у OMP).
Наблюдается та же закономерность: с ростом числа потоков эффективность снижается из-за ограниченной пропускной 


## 9. Выводы

- Версия TBB обеспечила максимальное ускорение **2.91×** (на 16 потоках), немного превзойдя OMP (2.67×) благодаря 
автоматической балансировке нагрузки.
- Наилучший баланс между производительностью и эффективностью достигается при 4 потоках: ускорение 2.74×, 
эффективность 68.5%.
- TBB удобен тем, что не требует ручной настройки grainsize или планировщика — `auto_partitioner` успешно справляется 
с распределением работы.
- В сравнении с OMP, TBB показал чуть более ровное поведение без резких просадок.
