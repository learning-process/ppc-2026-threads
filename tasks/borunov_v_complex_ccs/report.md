# Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – столбцовый (CCS).

- Студент: Борунов Владислав Алексеевич
- Группа: 3823Б1ПР3
- Вариант: 7
- Задача: Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – столбцовый (CCS).
- Технологии: SEQ, OMP, TBB, STL, ALL (OMP + TBB + STL в одном процессе)

Подробные отчёты по реализациям: [seq/report.md](seq/report.md), [omp/report.md](omp/report.md), [tbb/report.md](tbb/report.md), [stl/report.md](stl/report.md), [all/report.md](all/report.md).

## 1. Введение

Умножение разреженных матриц \(C = A \cdot B\) при больших размерностях и малой плотности эффективнее выполнять в сжатом столбцовом формате (CCS - Compressed Column Storage), не материализуя плотные массивы. Для комплексных коэффициентов каждая ненулевая запись хранится как std::complex<double>.

Цель работы: реализовать корректный последовательный baseline и параллельные версии (OpenMP, TBB, std::thread, комбинированный ALL), проверить их тестами фреймворка PPC и сравнить время выполнения на одном наборе perf-данных.

## 2. Постановка задачи

Вход: пара матриц (A, B) в CCS:

| Поле                      | Смысл                           |
|---------------------------|---------------------------------|
| num_rows, num_cols        | Размерность                     |
| col_ptrs[j]-col_ptrs[j+1] | Диапазон ненулевых в столбце j  |
| row_indices[k], values[k] | Строка и значение k-й ненулевой |

Выход: одна матрица C в CCS.

Ограничения: A.num_cols == B.num_rows; размер col_ptrs равен num_cols + 1.

Алгоритм (общий для всех реализаций): для каждого столбца j матрицы B накапливается столбец C[:,j] через произведения столбцов A[:,p] на скаляры B[p,j]. Ненулевые строки отслеживаются маркером или массивом флагов; малые по модулю значения отбрасываются (порог 1e-9).

## 3. Структура проекта

```
borunov_v_complex_ccs/
  common/include/common.hpp - SparseMatrix, InType, OutType
  seq/   — последовательная реализация
  omp/   — OpenMP
  tbb/   — Intel oneTBB
  stl/   — std::thread
  all/   — OMP + TBB + STL (гибрид)
  tests/functional/main.cpp
  tests/performance/main.cpp
  report.md, seq/report.md, …
```

## 4. Схема распараллеливания (обзор)

| Технология | Единица декомпозиции                         | Синхронизация                                                    |
|------------|----------------------------------------------|------------------------------------------------------------------|
| SEQ        | -                                            | нет                                                              |
| OMP        | столбцы B                                    | #pragma omp parallel, локальные буферы, последовательный merge   |
| TBB        | слоты tid (диапазоны столбцов)               | task_arena, parallel_for, merge вне параллельной фазы            |
| STL        | столбцы B                                    | std::thread, join() до merge                                     |
| ALL        | первая половина столбцов - OMP, вторая - TBB | барьер OMP, завершение task_arena, последовательная сборка (STL) |

## 5. Метрики (единые определения)

Во всех отчётах и таблицах:

| Метрика    | Определение                                                                                                                      |
|------------|----------------------------------------------------------------------------------------------------------------------------------|
| time       | Время этапа perf-теста, с: task_run - только RunImpl; pipeline - полный конвейер задачи                                          |
| speedup    | \(S = T_{\mathrm{seq,task\_run}} / T_{\mathrm{mode}}\); baseline - SEQ task_run                                                  |
| efficiency | \(E = S / \mathrm{workers}\) (в таблицах - в процентах)                                                                          |
| workers    | Число потоков, участвующих в RunImpl: OMP/TBB/ALL - PPC_NUM_THREADS (GetNumThreads()); STL - std::thread::hardware_concurrency() |

## 6. Экспериментальная установка

- CPU: i7-12650H
- RAM: 16 ГБ DDR5
- OS: Windows 11 
- Компилятор: MSVC 19.50 (Release сборка)

- Данные perf: m = k = n = 20000, в каждом столбце ~20 ненулевых, шаг по строкам max(1, num_rows/20) (tests/performance/main.cpp).
- Функциональные тесты: размеры (10,10,10), (20,15,25), (5,30,5), случайная разреженность 0.2, эталон - плотное умножение (tests/functional/main.cpp).
- Таймер: std::chrono::high_resolution_clock (настройка в BorunovVRunPerfTestThreads::SetPerfAttributes).

## 7. Результаты производительности

Baseline: SEQ task_run = 0.0408020800 с., SEQ pipeline = 0.0332925600 с.

### 7.1 task_run (только RunImpl)

Time:
| Workers | SEQ    | OMP    | TBB    | STL    | ALL    |
|---------|--------|--------|--------|--------|--------|
| 1       | 0.0408 | -      | -      | -      | -      |
| 4       | -      | 0.0618 | 0.0390 | 0.0580 | 0.0522 |
| 8       | -      | 0.0530 | 0.0474 | 0.0592 | 0.0425 |
| 16      | -      | 0.0437 | 0.0371 | 0.0494 | 0.0551 |

Speedup:
| Workers | OMP  | TBB  | STL  | ALL  |
|---------|------|------|------|------|
| 4       | 0.66 | 1.05 | 0.70 | 0.78 |
| 8       | 0.77 | 0.86 | 0.69 | 0.96 |
| 16      | 0.93 | 1.10 | 0.83 | 0.74 |

Efficiency:
| Workers | OMP   | TBB   | STL   | ALL   |
|---------|-------|-------|-------|-------|
| 4       | 16.5% | 26.2% | 17.6% | 19.5% |
| 8       | 9.6%  | 10.8% | 8.6%  | 12.0% |
| 16      | 5.8%  | 6.9%  | 5.2%  | 4.6%  |

### 7.2 pipeline (полный конвейер)

Time:
| Workers | SEQ    | OMP    | TBB    | STL    | ALL    |
|---------|--------|--------|--------|--------|--------|
| 1       | 0.0333 | -      | -      | -      | -      |
| 4       | -      | 0.0485 | 0.0513 | 0.0353 | 0.0356 |
| 8       | -      | 0.0230 | 0.0334 | 0.0303 | 0.0240 |
| 16      | -      | 0.0204 | 0.0525 | 0.0334 | 0.0345 |

Speedup:
| Workers | OMP  | TBB  | STL  | ALL  |
|---------|------|------|------|------|
| 4       | 0.69 | 0.65 | 0.94 | 0.94 |
| 8       | 1.45 | 1.00 | 1.10 | 1.39 |
| 16      | 1.63 | 0.63 | 1.00 | 0.96 |

Efficiency:
| Workers | OMP   | TBB   | STL   | ALL   |
|---------|-------|-------|-------|-------|
| 4       | 17.2% | 16.2% | 23.6% | 23.4% |
| 8       | 18.1% | 12.4% | 13.8% | 17.4% |
| 16      | 10.2% | 4.0%  | 6.2%  | 6.0%  |

### 7.3 Анализ (по таблицам выше)

| Режим    | Технология | Workers | speedup | efficiency |
|----------|------------|---------|---------|------------|
| task_run | TBB        | 4       | 1.05    | 26.2%      |
| task_run | TBB        | 16      | 1.10    | 6.9%       |
| pipeline | OMP        | 8       | 1.45    | 18.1%      |
| pipeline | OMP        | 16      | 1.63    | 10.2%      |
| pipeline | ALL        | 8       | 1.39    | 17.4%      |
| pipeline | STL        | 4       | 0.94    | 23.6%      |

Наблюдения:

- TBB в режиме task_run даёт \(S \ge 1\) при 4 и 16 потоках.
- OMP в режиме pipeline даёт наибольшее ускорение при 8-16 потоках (\(S \approx 1.45\)-\(1.63\)).
- ALL близок к лучшим режимам OMP/TBB на 8 потоках (pipeline, \(S=1.39\)), но не превосходит лучший одиночный вариант (OMP \(S=1.45\)).
- С ростом workers efficiency падает - накладные расходы merge и память на поток.

## 8. Сборка и воспроизведение

Из корня репозитория ppc-2026-threads:

```
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
```

Функциональные тесты задачи:

```bash
./build/bin/ mpiexec -n 8 .\ppc_perf_tests.exe --gtest_filter="*borunov_v_complex_ccs*"
```

Производительность (пример для 4, 8, 16 потоков):

```bash
./build/bin/ mpiexec -n 4 .\ppc_perf_tests.exe --gtest_filter="*borunov_v_complex_ccs*"
./build/bin/ mpiexec -n 8 .\ppc_perf_tests.exe --gtest_filter="*borunov_v_complex_ccs*"
./build/bin/ mpiexec -n 16 .\ppc_perf_tests.exe --gtest_filter="*borunov_v_complex_ccs*"
```


## 9. Корректность

Набор BorunovVRunFuncTestsThreads (tests/functional/main.cpp) сравнивает результат с плотным эталоном (допуск по значениям 1e-6). Локально пройдены реализации SEQ, OMP, TBB, STL, ALL (15 тестов: 3 размера x 5 технологий).

## 10. Выводы

Реализованы пять вариантов умножения CCS для комплексных матриц. Параллелизм - по столбцам B с локальными аккумуляторами и отдельной фазой слияния CCS.

На данных 20000x20000 (~20 ненулевых на столбец):

- TBB - лучший task_run (\(S \approx 1.05\)-\(1.10\) при 4-16 workers).
- OMP - лучший pipeline при 8 workers (\(S \approx 1.45\)).
- STL - устойчивый pipeline (\(S \approx 0.94\)-\(1.10\)).
- ALL - компромисс между OMP и TBB; на 8 workers pipeline близок к OMP, но чуть ниже лучшего одиночного результата.

Снижение efficiency с ростом числа потоков связано с последовательным merge и объёмом приватных структур \(O(\texttt{num\_rows})\) на поток.

## 11. Источники

1. Курс PPC, материалы по OpenMP, TBB — <https://learning-process.github.io/parallel_programming_slides/>
2. OpenMP Specification — <https://www.openmp.org/specifications/>
3. oneAPI TBB Documentation — <https://www.intel.com/content/www/us/en/docs/onetbb/developer-guide/current/overview.html>
