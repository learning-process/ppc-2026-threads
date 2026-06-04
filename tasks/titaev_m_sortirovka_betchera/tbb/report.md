# Поразрядная сортировка для вещественных чисел с чётно-нечётным слиянием Бэтчера - TBB

- Student: Титаев М., group 3823Б1ФИ1
- Technology: TBB
- Variant: 20

## 1. Introduction

TBB-версия использует Intel oneAPI Threading Building Blocks для параллельного выполнения этапов
преобразования данных и сети чётно-нечётного слияния Бэтчера. Цель — оценить производительность
библиотеки задач TBB на данной постановке и сравнить её с baseline.

## 2. Problem Statement

Задача совпадает с остальными версиями: отсортировать массив `InType = std::vector<double>` по
неубыванию, результат имеет тип `OutType = std::vector<double>`.

## 3. Baseline Algorithm (Sequential)

Baseline описан в `seq/report.md`. Он выполняет преобразование `double` в упорядоченные 64-битные ключи,
поразрядную сортировку по основанию 256 и обратное преобразование, после чего при размере, равном степени
двойки, применяет сеть Бэтчера.

## 4. Parallelization Scheme

TBB-реализация распараллеливает следующие этапы:

1. **Преобразование ключей:** прямое и обратное преобразование выполняется через
`oneapi::tbb::parallel_for` по `blocked_range<size_t>`.
2. **Поразрядная сортировка:** проходы счётной сортировки выполняются последовательно, что сохраняет
стабильность распределения.
3. **Сеть Бэтчера:** каждый шаг сравнения распараллелен через `oneapi::tbb::parallel_for`; пары индексов
на одном шаге не пересекаются.

Число потоков задаётся через `PPC_NUM_THREADS`.

## 5. Implementation Details

- Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`.
- Класс: `TitaevSortirovkaBetcheraTBB`.
- Параллельные участки выражены через `parallel_for` и `blocked_range`.
- Проходы поразрядной сортировки оставлены последовательными ради стабильности.

Память: один временный массив ключей и массив счётчиков на 256 корзин.

## 6. Experimental Setup

**Аппаратное обеспечение:**

- **CPU:** 12th Gen Intel(R) Core(TM) i5-12450H (2.00 GHz, 8 ядер / 12 потоков)
- **RAM:** 16 ГБ
- **OS:** Windows 11 Pro x64
- **MPI:** Microsoft MPI (MS-MPI) 10.1

**Инструменты:**

- **Сборка:** CMake
- **Компилятор:** MSVC 19.x
- **Конфигурация:** Release

**Окружение:**

- **PPC_NUM_THREADS:** задаёт число потоков для OMP, TBB, STL и потоковой части ALL.
- **PPC_NUM_PROC / mpiexec -n:** задаёт число MPI-процессов для ALL.
- Для ALL конфигурация записывается в формате `ranks × threads`.

**Генерация данных:**

- Тесты генерируют входные данные автоматически.
- Для performance-теста используется массив из `2^20` элементов с убывающими значениями.
- Внешние файлы с данными не используются.

## 7. Results and Discussion

### 7.1 Correctness

Корректность проверялась функциональными тестами из `tests/functional/main.cpp`. Выходной массив
проверяется на упорядоченность по неубыванию. Результаты TBB-версии совпадают с baseline на всех тестовых
наборах.

### 7.2 Performance

Используемые обозначения:

```txt
time — время выполнения performance-теста;
speedup = time_seq / time_mode;
efficiency = speedup / workers;
workers — количество исполнителей: потоков для OMP/TBB/STL, ranks × threads для ALL.
```

| Mode | Count | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- |
| seq | 1 | 0.159970 | 1.00 | N/A |
| tbb | 2 | 0.123646 | 1.29 | 64.69% |
| tbb | 4 | 0.074043 | 2.16 | 54.01% |
| tbb | 8 | 0.063984 | 2.50 | 31.25% |

Производительность TBB растёт с числом потоков и приближается к результату OMP при восьми потоках.
Накладные расходы на планирование задач TBB заметны на малом числе потоков, поэтому на двух потоках
ускорение невелико.

## 8. Conclusions

TBB-версия показывает рост производительности с увеличением числа потоков, наилучший результат получен на
восьми потоках. Эффективность снижается с ростом числа потоков из-за накладных расходов на планирование
задач и большого числа синхронизированных проходов в сети Бэтчера.

## 9. References

1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.
2. oneAPI Threading Building Blocks Documentation.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: `std::thread`.

## Appendix (Optional)

```cpp
void TitaevSortirovkaBetcheraTBB::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
      oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, n),
                                [&](const oneapi::tbb::blocked_range<size_t> &r) {
                                  for (size_t i = r.begin(); i < r.end(); i++) {
                                    const size_t l = i ^ j;
                                    if (l > i && l < n) {
                                      const bool ascending = ((i & k) == 0);
                                      const bool need_swap =
                                          ascending ? (result[i] > result[l]) : (result[i] < result[l]);
                                      if (need_swap) {
                                        std::swap(result[i], result[l]);
                                      }
                                    }
                                  }
                                });
    }
  }
}
```
