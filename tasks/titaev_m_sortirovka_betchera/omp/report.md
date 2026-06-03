# Поразрядная сортировка для вещественных чисел с чётно-нечётным слиянием Бэтчера - OMP

- Student: Титаев М., group 3823Б1ФИ1
- Technology: OMP
- Variant: 20

## 1. Introduction

OMP-версия ускоряет поразрядную сортировку массива `double` и сеть чётно-нечётного слияния Бэтчера за
счёт многопоточности OpenMP. Цель — распараллелить наиболее затратные по времени этапы и сравнить
производительность с последовательным baseline.

## 2. Problem Statement

Задача совпадает с остальными версиями: отсортировать массив `InType = std::vector<double>` по
неубыванию, результат имеет тип `OutType = std::vector<double>`.

## 3. Baseline Algorithm (Sequential)

Baseline описан в `seq/report.md`. Он выполняет преобразование `double` в упорядоченные 64-битные ключи,
поразрядную сортировку по основанию 256 и обратное преобразование, после чего при размере, равном степени
двойки, применяет сеть Бэтчера.

## 4. Parallelization Scheme

OMP-реализация распараллеливает следующие этапы:

1. **Преобразование ключей:** цикл преобразования `double` в `uint64_t` распараллелен через
`#pragma omp parallel for`.
2. **Поразрядная сортировка:** на каждом проходе каждый поток заполняет собственную локальную гистограмму
корзин, затем гистограммы суммируются, и выполняется стабильное распределение.
3. **Сеть Бэтчера:** на каждом шаге сравнения пары индексов обрабатываются параллельным циклом
`#pragma omp parallel for`; пары на одном шаге не пересекаются, поэтому гонок нет.

Число потоков задаётся через `PPC_NUM_THREADS`.

## 5. Implementation Details

- Файлы: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`.
- Класс: `TitaevSortirovkaBetcheraOMP`.
- Локальные гистограммы исключают гонки при подсчёте корзин.
- Распределение элементов по корзинам выполняется последовательно для сохранения стабильности.

Память: один временный массив ключей и по одной гистограмме на каждый поток.

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
проверяется на упорядоченность по неубыванию. Результаты OMP-версии совпадают с baseline на всех тестовых
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
| omp | 2 | 0.056191 | 2.85 | 142.34% |
| omp | 4 | 0.063745 | 2.51 | 62.73% |
| omp | 8 | 0.080437 | 1.99 | 24.86% |

Лучший результат достигается на двух потоках. С ростом числа потоков время увеличивается: сеть Бэтчера
выполняет порядка нескольких сотен синхронизированных проходов по массиву, и накладные расходы на барьеры
между проходами начинают преобладать над выигрышем от параллелизма.

## 8. Conclusions

OMP-версия ускоряет сортировку относительно baseline, наибольшее ускорение получено на двух потоках.
Дальнейшее увеличение числа потоков снижает эффективность из-за большого количества синхронизаций в сети
Бэтчера.

## 9. References

1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.
2. oneAPI Threading Building Blocks Documentation.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: `std::thread`.

## Appendix (Optional)

```cpp
void TitaevSortirovkaBetcheraOMP::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
#pragma omp parallel for
      for (long long ii = 0; ii < static_cast<long long>(n); ii++) {
        const size_t i = static_cast<size_t>(ii);
        const size_t l = i ^ j;
        if (l > i) {
          const bool ascending = ((i & k) == 0);
          const bool need_swap = ascending ? (result[i] > result[l]) : (result[i] < result[l]);
          if (need_swap) {
            std::swap(result[i], result[l]);
          }
        }
      }
    }
  }
}
```
