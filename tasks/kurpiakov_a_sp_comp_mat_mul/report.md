# Умножение разреженных комплексных матриц в CSR

- Student: Курпяков Алексей Георгиевич, group 3823Б1ФИ3  
- Variant: 6  
- Local reports: seq/report.md, omp/report.md, tbb/report.md, stl/report.md, all/report.md

## 1. Введение

Отчет объединяет результаты по SEQ/OMP/TBB/STL/ALL на одном вычислительном ядре и одном наборе входных данных.

## 2. Единая постановка задачи

Вход — две CSR‑матрицы с комплексными элементами, выход — их произведение в CSR. Корректность оценивается сравнением с SEQ‑эталоном.

## 3. Единая методика эксперимента

Окружение: AMD Ryzen 7 7840HS, 4 GB RAM, Ubuntu 22.04, GCC 13, Release.  
Переменные: `PPC_NUM_THREADS = 4`, `PPC_NUM_PROC = 4`.  
Данные perf‑теста: трехдиагональные матрицы размера $N = 5{,}000{,}000$.

Кодовый фрагмент генерации данных (tasks/kurpiakov_a_sp_comp_mat_mul/tests/performance/main.cpp):

```cpp
constexpr int kSize = 5000000;

auto a = MakeTridiagonal(kSize);
auto b = MakeTridiagonal(kSize);

expected_output_ = a.Multiply(b);
```

**Определения метрик:**

- `time` — значение из perf‑логов (условные единицы, без пересчета).  
- $S = T_{seq}/T_{x}$ — ускорение для того же режима.  
- `workers` $W$ — число рабочих единиц (SEQ: 1, OMP/TBB/STL: число потоков, ALL: $ranks \times threads$).  
- $E = S / W$ — эффективность.

## 4. Сводка корректности

Функциональные тесты выполнялись локально для SEQ в tasks/kurpiakov_a_sp_comp_mat_mul/tests/functional/main.cpp.  
Для всех backend‑ов корректность проверялась локально в perf‑тестах сравнением с эталоном `Multiply`.

## 5. Агрегированные результаты

### pipeline

| Backend | Workers | Time | Speedup | Efficiency |
| --- | ---: | ---: | ---: | ---: |
| SEQ | 1 | 8125 | 1.000 | 1.000 |
| OMP | 4 | 9358 | 0.868 | 0.217 |
| TBB | 4 | 7438 | 1.092 | 0.273 |
| STL | 4 | 9701 | 0.838 | 0.210 |
| ALL | 16 | 13759 | 0.590 | 0.037 |

### task_run

| Backend | Workers | Time | Speedup | Efficiency |
| --- | ---: | ---: | ---: | ---: |
| SEQ | 1 | 8509 | 1.000 | 1.000 |
| OMP | 4 | 9036 | 0.942 | 0.236 |
| TBB | 4 | 6843 | 1.243 | 0.311 |
| STL | 4 | 8545 | 0.996 | 0.249 |
| ALL | 16 | 12695 | 0.670 | 0.042 |

## 6. Интерпретация различий

TBB показывает наименьшее время в обоих режимах; OMP и STL не ускоряют SEQ на
выбранном размере. Для ALL снижение производительности согласуется с
дополнительными коммуникациями, но без профилирования это остается гипотезой.

## 7. Репродуцируемость

Сборка:

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Функциональные тесты:

```bash
PPC_NUM_THREADS=4 ./build/bin/ppc_func_tests --gtest_filter='*kurpiakov*'
```

Perf‑тесты:

```bash
PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
 --gtest_filter='*kurpiakov*<тип задачи>*'
PPC_NUM_PROC=4 PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
 --gtest_filter='*kurpiakov*all*'
```

## 8. Графики

Графики не использовались.

## 9. Заключение

На выбранных входах и конфигурации наибольшее ускорение дает TBB; результаты приведены в таблицах выше.
