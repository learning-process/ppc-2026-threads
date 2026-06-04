# Умножение плотных матриц. Элементы типа `double`. Алгоритм Штрассена

- Студент: Ахметов Даниил Данисович, группа 3823Б1ПР2
- Технология: OMP
- Вариант: 3

## 1. Введение

После `SEQ` задача перенесена на OpenMP. Для Штрассена это естественный
шаг: внутри алгоритма много независимых участков работы. В этой реализации
OpenMP применяется внутри итеративного стека, а не только на верхнем
уровне `M1..M7`.

Общий контекст - в [report.md](../report.md), baseline - в
[seq/report.md](../seq/report.md).

## 2. Постановка задачи

Та же постановка, что и в `SEQ`: вычислить `C = A * B` для плотных
матриц `double`. Корректность - относительно наивного эталона
с допуском `1e-7`.

## 3. Описание алгоритма (базового/последовательного)

Математическая часть совпадает с `SEQ`: итеративный Штрассен,
`kThreshold = 64`, padding до степени двойки. Подробнее - в
[seq/report.md](../seq/report.md).

## 4. Схема распараллеливания

**Декомпозиция:** параллелизация внутри итеративного стека.

**Где используется OpenMP:**

- `StandardMultiply` - базовое умножение;
- `Split`, `Merge` - разбиение и сборка подматриц;
- `AddInto`, `SubInto` - поэлементные операции;
- `PadTo` и копирование результата после padding.

**Планирование и атрибуты:**

- `#pragma omp parallel for`;
- `default(none)`, явное перечисление переменных;
- `schedule(static)`;
- `if (size >= kParallelThreshold)`, `kParallelThreshold = 256`.

**Синхронизация:**

- неявные барьеры в конце parallel for;
- `reduction`/`atomic`/`critical` не нужны - потоки пишут
  в разные строки/элементы.

## 5. Детали реализации

**Файлы:**

- `omp/include/ops_omp.hpp`
- `omp/src/ops_omp.cpp`

**Изменения относительно SEQ:**

- сохранён стек `Frame`;
- добавлены OpenMP-директивы в горячие участки;
- число потоков - `ppc::util::GetNumThreads()`.

**Пайплайн:** `ValidationImpl` -> `PreProcessingImpl` -> `RunImpl` ->
`PostProcessingImpl`.

## 6. Экспериментальная установка

- CPU: AMD Ryzen 7 2700, 8 ядер / 16 потоков;
- RAM: 16 GB;
- OS: Windows 10 x64;
- IDE: Visual Studio Code;
- compiler: MSVC (сборка через CMake);
- build type: `Release`;
- `PPC_NUM_THREADS=4`, `OMP_NUM_THREADS=4`.

**Данные:** performance - `1024 × 1024`, случайные `[-10, 10]`.

**Команды:**

```powershell
cd build\bin
$env:OMP_NUM_THREADS = "4"
$env:PPC_NUM_THREADS = "4"
.\ppc_func_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*_omp_*
.\ppc_perf_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double_omp_enabled*
```

## 7. Результаты и обсуждение

### 7.1 Корректность

- размеры `64`, `128`, `256`;
- сравнение с наивным эталоном, допуск `1e-7`;
- локально все функциональные тесты пройдены.

### 7.2 Производительность

| Режим | Число потоков | Время, с | Ускорение | Эффективность |
| --- | ---: | ---: | ---: | ---: |
| seq | 1 | 0.938515 | 1.00 | N/A |
| omp | 4 | 1.705478 | 0.55 | 13.8% |

На `1024 × 1024` `OMP` оказался медленнее `SEQ`. Вероятная причина -
накладные расходы OpenMP на множество относительно небольших циклов
в итеративной реализации.

## 8. Заключение

`OMP`-версия корректна и используется как основа для `ALL`, но на
текущей конфигурации не показала ускорения. Для этой задачи эффективнее
оказалась параллелизация верхнего уровня `M1..M7` (`TBB`, `STL`).

## 9. Источники

1. [OpenMP Specifications](https://www.openmp.org/specifications/)
2. [Parallel Programming Course](https://learning-process.github.io/parallel_programming_course/ru/common_information/threading_tasks.html)

## Приложение (опционально)

```cpp
#pragma omp parallel for default(none) shared(a, b, c, size) \
    schedule(static) if (size >= kParallelThreshold)
for (size_t i = 0; i < size; ++i) { /* ... */ }
```
