# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) - STL

- Student: Тельнов А., group 3823Б1ФИ1
- Technology: STL
- Variant: 9

## 1. Introduction

STL-версия использует `std::thread` для ручного управления потоками. Эта реализация показывает, как метод прямоугольников можно распараллелить без OpenMP и TBB, явно разделив диапазон индексов между потоками.

## 2. Problem Statement

Нужно вычислить интеграл функции `f(x1, ..., xd) = x1 + ... + xd` на области `[0, 1]^d`. Вход — `n` и `dimensions`, выход — значение интеграла типа `double`.

## 3. Baseline Algorithm (Sequential)

Baseline описан в `seq/report.md`. Последовательный алгоритм обходит все точки сетки в одном потоке и используется для расчёта ускорения.

## 4. Parallelization Scheme

Общее число точек `total_points = n^dimensions` делится на `thread_count` непересекающихся блоков. Размеры блоков отличаются не более чем на один элемент:

```cpp
block = total_points / thread_count;
remainder = total_points % thread_count;
```

Каждый поток вычисляет сумму на своём диапазоне `[begin, end)` и записывает результат в `partial_sums[i]`.

Важный момент: сначала запускаются все потоки в цикле:

```cpp
for (int i = 0; i < thread_count; ++i) {
  threads[i] = std::thread(...);
}
```

И только после запуска всех потоков выполняется ожидание завершения:

```cpp
for (auto &thread : threads) {
  thread.join();
}
```

Такой порядок обеспечивает реальный параллелизм. Если бы `join()` вызывался сразу после создания каждого потока, выполнение стало бы почти последовательным.

Синхронизация выполняется через `join()`. Мьютексы не используются, потому что каждый поток записывает только в свою ячейку `partial_sums[i]`.

## 5. Implementation Details

- Файлы: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`.
- Класс: `TelnovAIntegralRectangleSTL`.
- Число потоков берётся через `ppc::util::GetNumThreads()` и ограничивается числом точек.
- Лямбда `calculate_part` вычисляет частичную сумму на заданном диапазоне.
- После `join()` главный поток последовательно складывает значения из `partial_sums` и умножает результат на `h^dimensions`.

Память: `std::vector<std::thread>` и `std::vector<double> partial_sums`. Размер дополнительной памяти пропорционален числу потоков, а не числу точек сетки.

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
- Для performance-теста используется вход `InType{50, 4}`: 50 разбиений по каждой координате и размерность 4.
- Внешние файлы с данными не используются.

## 7. Results and Discussion

### 7.1 Correctness

Корректность проверялась функциональными тестами из `tests/functional/main.cpp`. Для набора параметров `(n, dimensions)` вычисленное значение сравнивается с аналитическим результатом:

```txt
I = dimensions / 2
```

Это значение получается для функции `f(x1, ..., xd) = x1 + ... + xd` на единичном гиперкубе `[0, 1]^d`: интеграл каждой координаты равен `1/2`, поэтому сумма по `d` координатам равна `d / 2`. В тестах используется допустимое отклонение, зависящее от числа разбиений `n`.

### 7.2 Performance

Используемые обозначения:

```txt
time — время выполнения performance-теста;
speedup = time_seq / time_mode;
efficiency = speedup / workers;
workers — количество исполнителей: потоков для OMP/TBB/STL, ranks × threads для ALL.
```
| Mode | Count | Time, s | Speedup | Efficiency |
|-------------|-------|---------|---------|------------|
| seq | 1 | 0.187439 | 1.00 | N/A |
| stl | 2 | 0.098049 | 1.91 | 95.58% |
| stl | 4 | 0.070281 | 2.67 | 66.68% |
| stl | 8 | 0.054564 | 3.43 | 42.94% |

## 8. Conclusions

Реализация вычисляет многомерный интеграл методом прямоугольников без хранения полной сетки в памяти. Последовательная версия используется как baseline, а параллельные версии сравниваются с ней по времени выполнения, ускорению и эффективности. Основное ограничение масштабируемости связано с ростом числа точек `n^dimensions`, накладными расходами на управление параллелизмом и редукцией частичных сумм.

## 9. References

1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.
2. oneAPI Threading Building Blocks Documentation.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: `std::thread`.

## Appendix (Optional)

Основной фрагмент `RunImpl()`:

```cpp
bool TelnovAIntegralRectangleSTL::RunImpl() {
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / static_cast<double>(n);

  const auto total_points = static_cast<int64_t>(std::pow(n, d));

  int thread_count = ppc::util::GetNumThreads();
  thread_count = std::max(1, std::min(thread_count, static_cast<int>(total_points)));

  std::vector<std::thread> threads(thread_count);
  std::vector<double> partial_sums(thread_count, 0.0);

  const int64_t block = total_points / thread_count;
  const int64_t remainder = total_points % thread_count;

  auto calculate_part = [n, d, a, h](int64_t begin, int64_t end) {
    double local_sum = 0.0;

    for (int64_t idx = begin; idx < end; ++idx) {
      int64_t current = idx;
      double f_value = 0.0;

      for (int dim = 0; dim < d; ++dim) {
        const int coord_index = static_cast<int>(current % n);
        current /= n;

        const double x = a + ((static_cast<double>(coord_index) + 0.5) * h);
        f_value += x;
      }

      local_sum += f_value;
    }

    return local_sum;
  };

  int64_t begin = 0;
  for (int i = 0; i < thread_count; ++i) {
    const int64_t current_block = block + (i < remainder ? 1 : 0);
    const int64_t end = begin + current_block;

    threads[i] = std::thread([&, i, begin, end]() { partial_sums[i] = calculate_part(begin, end); });

    begin = end;
  }

  for (auto &thread : threads) {
    thread.join();
  }

  double result = 0.0;
  for (const auto &value : partial_sums) {
    result += value;
  }

  GetOutput() = result * std::pow(h, d);
  return true;
}
```
