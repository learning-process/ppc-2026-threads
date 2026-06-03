# Умножение плотных матриц. Элементы типа `double`. Алгоритм Штрассена

- Студент: Ахметов Даниил Данисович, группа 3823Б1ПР2
- Технология: STL
- Вариант: 3

## 1. Введение

Эта версия показывает распараллеливание на стандартных средствах C++
без OpenMP и без отдельного runtime TBB. Используются `std::async`
и `std::future` (backend называется `STL`, хотя в коде не `std::thread`
на верхнем уровне, а `std::async`).

Общий контекст - в [report.md](../report.md).

## 2. Постановка задачи

Та же постановка, что и в `SEQ`: `C = A * B`, допуск `1e-7`.
Подробнее - в [seq/report.md](../seq/report.md).

## 3. Описание алгоритма (базового/последовательного)

Внутренняя часть совпадает с `TBB`: `StrassenSeqImpl`, `kCutoff = 256`,
блочное умножение при малых размерах.

## 4. Схема распараллеливания

**Декомпозиция:** верхний уровень `M1..M7`.

**Примитив:**

- `std::async(std::launch::async, ...)` для каждой ветви;
- синхронизация через `future.get()` на главном потоке.

**Данные:**

- каждая задача пишет в свой буфер `m1..m7`;
- `mutex` и `atomic` в горячем участке не нужны.

**Планирование:** семь независимых асинхронных задач, по одной на
каждое промежуточное произведение.

## 5. Детали реализации

**Файлы:**

- `stl/include/ops_stl.hpp`
- `stl/src/ops_stl.cpp`

**Ключевые функции:**

- `StrassenTopStl(...)` - верхний уровень с `std::async`;
- `StrassenSeqImpl(...)` - внутренняя последовательная часть.

Матрицы хранятся в `std::vector<double>`.

## 6. Экспериментальная установка

- CPU: AMD Ryzen 7 2700, 8 ядер / 16 потоков;
- RAM: 16 GB;
- OS: Windows 10 x64;
- IDE: Visual Studio Code;
- compiler: MSVC (сборка через CMake);
- build type: `Release`;
- `PPC_NUM_THREADS=4`.

**Команды:**

```powershell
cd build\bin
$env:PPC_NUM_THREADS = "4"
.\ppc_func_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*_stl_*
.\ppc_perf_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double_stl_enabled*
```

## 7. Результаты и обсуждение

### 7.1 Корректность

- размеры `64`, `128`, `256`;
- сравнение с наивным эталоном;
- локально все тесты пройдены.

### 7.2 Производительность

| Режим | Число потоков | Время, с | Ускорение | Эффективность |
| --- | ---: | ---: | ---: | ---: |
| seq | 1 | 0.938515 | 1.00 | N/A |
| stl | 4 | 0.100045 | 9.38 | 234.5% |

`STL` показал лучшее локальное время и оказался очень близок к `TBB`.

## 8. Заключение

Верхний уровень Штрассена можно эффективно распараллелить
стандартными средствами C++. `STL` - лучший результат
по локальному времени среди всех backend-ов.

## 9. Источники

1. [cppreference: std::async](https://en.cppreference.com/w/cpp/thread/async)
2. [Parallel Programming Course](https://learning-process.github.io/parallel_programming_course/ru/common_information/threading_tasks.html)

## Приложение (опционально)

```cpp
auto f1 = std::async(std::launch::async, [&] { /* M1 */ });
// ... f2 .. f7 ...
f1.get(); /* ... */
```
