# Поразрядная сортировка для вещественных чисел с чётно-нечётным слиянием Бэтчера - STL

- Student: Титаев М., group 3823Б1ФИ1
- Technology: STL
- Variant: 20

## 1. Introduction

STL-версия реализует параллелизм с помощью стандартной библиотеки потоков `std::thread`. Цель — вручную
разбить работу между потоками на этапах преобразования данных и сети чётно-нечётного слияния Бэтчера и
сравнить производительность с baseline.

## 2. Problem Statement

Задача совпадает с остальными версиями: отсортировать массив `InType = std::vector<double>` по
неубыванию, результат имеет тип `OutType = std::vector<double>`.

## 3. Baseline Algorithm (Sequential)

Baseline описан в `seq/report.md`. Он выполняет преобразование `double` в упорядоченные 64-битные ключи,
поразрядную сортировку по основанию 256 и обратное преобразование, после чего при размере, равном степени
двойки, применяет сеть Бэтчера.

## 4. Parallelization Scheme

STL-реализация распараллеливает следующие этапы:

1. **Преобразование ключей:** диапазон индексов делится на непрерывные блоки, каждый обрабатывается
отдельным `std::thread`.
2. **Поразрядная сортировка:** на каждом проходе каждый поток заполняет собственную локальную гистограмму
корзин на своём блоке; гистограммы суммируются, после чего выполняется стабильное распределение.
3. **Сеть Бэтчера:** на каждом шаге сравнения диапазон индексов делится между потоками; пары индексов на
одном шаге не пересекаются, поэтому гонок нет.

Объединение результатов выполняется после `join()` всех потоков.

## 5. Implementation Details

- Файлы: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`.
- Класс: `TitaevSortirovkaBetcheraSTL`.
- Число потоков определяется через `std::thread::hardware_concurrency()`.
- Локальные гистограммы и поблочное разбиение исключают гонки данных.

Память: один временный массив ключей, по одной гистограмме на каждый поток и вектор объектов потоков.

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
проверяется на упорядоченность по неубыванию. Результаты STL-версии совпадают с baseline на всех тестовых
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
| stl | 2 | 0.179218 | 0.89 | 44.62% |
| stl | 4 | 0.164399 | 0.97 | 24.33% |
| stl | 8 | 0.165794 | 0.96 | 12.06% |

Время STL-версии слабо зависит от заданного числа потоков, поскольку число рабочих потоков определяется
аппаратной конфигурацией процессора, а не переменной окружения. На больших массивах сеть Бэтчера создаёт
потоки в каждом из множества проходов, и затраты на создание и объединение потоков превышают выигрыш от
распараллеливания, поэтому ускорение относительно baseline отсутствует.

## 8. Conclusions

STL-версия на данной постановке не даёт ускорения относительно baseline. Основная причина —
многократное создание потоков на каждом шаге сети Бэтчера и независимость числа потоков от заданного
параметра. Для улучшения масштабируемости потребовался бы пул потоков и сокращение числа точек
синхронизации.

## 9. References

1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.
2. oneAPI Threading Building Blocks Documentation.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: `std::thread`.

## Appendix (Optional)

```cpp
void TitaevSortirovkaBetcheraSTL::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  const unsigned int num_threads = GetThreadCount();

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
      std::vector<std::thread> threads;
      threads.reserve(num_threads);
      const size_t chunk = (n + num_threads - 1) / num_threads;

      for (unsigned int t = 0; t < num_threads; t++) {
        const size_t begin = t * chunk;
        const size_t end = std::min(begin + chunk, n);
        if (begin >= end) {
          break;
        }
        threads.emplace_back([&result, n, k, j, begin, end]() {
          for (size_t i = begin; i < end; i++) {
            const size_t l = i ^ j;
            if (l > i && l < n) {
              const bool ascending = ((i & k) == 0);
              const bool need_swap = ascending ? (result[i] > result[l]) : (result[i] < result[l]);
              if (need_swap) {
                std::swap(result[i], result[l]);
              }
            }
          }
        });
      }
      for (auto &th : threads) {
        th.join();
      }
    }
  }
}
```
