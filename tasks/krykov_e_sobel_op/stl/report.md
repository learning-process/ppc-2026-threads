# Выделение рёбер на изображении с использованием оператора Собеля - STL

- Student: Крюков Е.
- Technology: STL 
- Variant: 27

## 1. Контекст

`std::thread` - инструмент ручного управления потоками из стандартной библиотеки
C++11. В отличие от OMP и TBB, здесь программист явно создаёт потоки,
распределяет данные, определяет критические секции и вызывает `join`.
Это требует больше кода, но даёт полный контроль над декомпозицией.
Baseline - `seq/report.md`.

## 2. Постановка задачи

Идентична SEQ-версии: `Image{width, height, data}` - `std::vector<int>`.
Полное совпадение с SEQ на всём тестовом наборе.

## 3. Базовый алгоритм

Подробно описан в `seq/report.md`. RGB - grayscale - свёртка ядром Собеля 3×3
по внутренним пикселям - `magnitude = sqrt(Gx² + Gy²)`.

## 4. Схема распараллеливания

Вспомогательная функция `ComputeMagnitude` вынесена в анонимное пространство
имён - она вычисляет магнитуду градиента для одного пикселя.

```cpp
// File: stl/src/ops_stl.cpp - RunImpl (ключевые фрагменты)
unsigned int num_threads = std::thread::hardware_concurrency();
if (num_threads == 0) num_threads = 2;

const unsigned int total_rows = static_cast<unsigned int>(h - 2);
const unsigned int rows_per_thread = total_rows / num_threads;
const unsigned int remainder = total_rows % num_threads;

std::vector<std::thread> threads;
threads.reserve(num_threads);

unsigned int next_start = 1;
for (unsigned int i = 0; i < num_threads; ++i) {
  unsigned int chunk = rows_per_thread + (i < remainder ? 1 : 0);
  unsigned int start = next_start;
  unsigned int end = start + chunk;
  next_start = end;
  threads.emplace_back([&, start, end]() {
    for (unsigned int row = start; row < end; ++row)
      for (int col = 1; col < w - 1; ++col)
        output[(row * w) + col] =
            ComputeMagnitude(gray, w, (int)row, col, gx_kernel, gy_kernel);
  });
}

for (auto& t : threads) t.join();
```

**Разбиение диапазона:** внутренние строки `[1, h-2]` делятся между потоками
равными блоками; «лишние» строки распределяются по одной на первые `remainder`
потоков.

**Что делает каждый поток:** обрабатывает свой непрерывный диапазон строк,
записывает результаты в `output`. Диапазоны не пересекаются - гонок на запись
нет.

**`join`:** вызывается в отдельном цикле **после** запуска всех потоков.
Это принципиально: `join` внутри цикла создания потоков сериализовал бы их
выполнение.

**Синхронизация:** явных `mutex`, `atomic` или `condition_variable` не
используется - достаточно диапазонной независимости. Барьер обеспечивается
`join` для каждого потока перед возвратом из `RunImpl`.

## 5. Детали реализации

Файлы: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`.

Вспомогательная функция `ComputeMagnitude` вынесена в анонимный namespace:

```cpp
// File: stl/src/ops_stl.cpp - вспомогательная функция
int ComputeMagnitude(const std::vector<int>& gray, int w, int row, int col,
                     const std::array<std::array<int, 3>, 3>& gx_kernel,
                     const std::array<std::array<int, 3>, 3>& gy_kernel) {
  int gx = 0, gy = 0;
  for (int ky = -1; ky <= 1; ++ky)
    for (int kx = -1; kx <= 1; ++kx) {
      int pixel = gray[(row + ky) * w + (col + kx)];
      gx += pixel * gx_kernel.at(ky + 1).at(kx + 1);
      gy += pixel * gy_kernel.at(ky + 1).at(kx + 1);
    }
  return static_cast<int>(std::sqrt((double)(gx * gx + gy * gy)));
}
```

Передача ядер по const-ссылке безопасна: ядра объявлены на стеке `RunImpl`
до создания потоков и живут на всё время их работы.

`ValidationImpl`, `PreProcessingImpl`, `PostProcessingImpl` идентичны SEQ.

## 6. Проверка корректности

Все три функциональных теста пройдены. При `num_threads = 1` результат совпадает
с SEQ; при `num_threads = 4` и `num_threads = 8` результат не меняется.
Гонок на запись нет по конструкции (непересекающиеся диапазоны).
Чтение `gray` из нескольких потоков безопасно - объект не изменяется.

## 7. Экспериментальная среда

| Параметр         | Значение                             |
| ---------------- | ------------------------------------ |
| CPU              | AMD Ryzen 3 3100 (4 физических ядра) |
| ОС               | Windows 10                           |
| Компилятор       | MSVC (Release)                       |
| CMake build type | Release                              |
| Число потоков    | `hardware_concurrency()` (8 на тестовой машине) |

Команда запуска:

```bash
mpiexec -n 4 .\build\bin\ppc_perf_tests.exe --gtest_filter="*KrykovE*"
```

## 8. Результаты

| size | SEQ task_run, с | STL task_run, с | Speedup | Efficiency (8 th.) |
| ---- | --------------- | --------------- | ------- | ------------------ |
| 512  | 0.00738         | 0.00627         | 1.18x   | 14.7%              |
| 1024 | 0.03471         | 0.02205         | 1.57x   | 19.7%              |
| 2048 | 0.11435         | 0.07197         | 1.59x   | 19.9%              |

| size | SEQ pipeline, с | STL pipeline, с | Speedup |
| ---- | --------------- | --------------- | ------- |
| 512  | 0.00851         | 0.01088         | 0.78x   |
| 1024 | 0.03374         | 0.02032         | 1.66x   |
| 2048 | 0.12496         | 0.08187         | 1.53x   |

На размере 512 в `pipeline`-режиме STL проигрывает SEQ: overhead создания
8 потоков превышает выигрыш от параллелизации на малой задаче.
Начиная с 1024 ускорение становится ощутимым (~1.57-1.59x).

## 9. Выводы

Ручная реализация `std::thread` показала ускорение, сопоставимое с OMP при
больших размерах (1024, 2048), но более высокие накладные расходы на малых
задачах. Ключевое достоинство - полный контроль над декомпозицией; недостаток -
существенно больше кода и ответственность программиста за корректность разбиения
и порядок `join`.
