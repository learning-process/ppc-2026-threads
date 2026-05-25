# Умножение плотных матриц. Элементы типа double. Блочная схема, алгоритм Кэннона. — TBB

- Student: Черемхин Андрей Александрович
- Technology: TBB
- Variant: 1

## 1. Контекст

oneTBB-версия использует ту же блочную декомпозицию, что и SEQ/OMP, но планирование независимых работ
передаётся runtime-библиотеке TBB. Это позволяет описать работу как набор диапазонов, а не вручную управлять
потоками.

## 2. Постановка задачи

Требуется вычислить произведение двух квадратных матриц `C = A * B`. Вход и выход совпадают с общей
постановкой: `n`, вектор `A`, вектор `B`, результат в построчном одномерном виде.

## 3. Базовый алгоритм

Алгоритм выбирает `q`, `bs` и `np`, копирует данные в padded-буферы, вычисляет блоки результата и копирует
ответ обратно. Независимость блоков `C[bi, bj]` позволяет использовать двумерный диапазон TBB.

## 4. Схема распараллеливания

Конкуренция ограничивается объектом:

```cpp
oneapi::tbb::global_control control(oneapi::tbb::global_control::max_allowed_parallelism, requested_threads);
```

Здесь `requested_threads` берётся из `PPC_NUM_THREADS`. Это ограничивает максимальный параллелизм TBB runtime
на время выполнения `RunImpl`.

Для копирования входа и выхода используется `oneapi::tbb::parallel_for` по строкам. Основное вычисление
использует `oneapi::tbb::blocked_range2d<std::int64_t>(0, q64, 0, q64)`:

```cpp
oneapi::tbb::parallel_for(oneapi::tbb::blocked_range2d<std::int64_t>(0, q64, 0, q64),
                          [&](const oneapi::tbb::blocked_range2d<std::int64_t> &range) {
  for (std::int64_t bi = range.rows().begin(); bi != range.rows().end(); ++bi) {
    for (std::int64_t bj = range.cols().begin(); bj != range.cols().end(); ++bj) {
      for (std::size_t step = 0; step < q; ++step) {
        const std::size_t bk = (static_cast<std::size_t>(bi) + static_cast<std::size_t>(bj) + step) % q;
        MulAddBlock(a, b, c, np, bs, static_cast<std::size_t>(bi), bk, static_cast<std::size_t>(bj));
      }
    }
  }
});
```

Grain size явно не задан, поэтому используется разбиение по умолчанию. Это упрощает код, но оставляет
runtime-у право выбирать размер подзадач. Для очень малых `q` накладные расходы TBB могут быть заметны.

## 5. Детали реализации

Файлы реализации: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`.

Локальных аккумуляторов не требуется: каждый TBB-task записывает в свой блок `C[bi, bj]`. Общие входные
матрицы `a` и `b` только читаются после фазы копирования. Запись в `out` также разделена по строкам, поэтому
гонки отсутствуют.

## 6. Проверка корректности

TBB-версия проверяется теми же функциональными тестами, что и остальные реализации. Результат сравнивается с
независимым последовательным oracle с допуском `1e-7`.

Команда для потоковых тестов:

```bash
cmake -S . -B build -D USE_COVERAGE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_BUILD_TYPE=Releas
cmake --build build -j --parallel
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
mpirun --oversubscribe -n 1 ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_tbb*'
mpirun --oversubscribe -n 2 ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_tbb*'
mpirun --oversubscribe -n 4 ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_tbb*'
```

## 7. Экспериментальная среда

- OS: Linux 6.6.114.1-microsoft-standard-WSL2 x86_64;
- CPU: AMD Ryzen 5 5600, 6 cores / 12 threads;
- Compiler: GCC 13.3.0;
- Build type: `Release`;
- TBB concurrency control: `global_control::max_allowed_parallelism`;
- Worker count source: `PPC_NUM_THREADS`.

## 8. Результаты

| size | workers | mode | median time, s | speedup vs seq | efficiency |
| ---: | ---: | --- | ---: | ---: | ---: |
| 640 | 1 | task | 0.8220091702 | 2.9644 | 2.9644 |
| 640 | 2 | task | 0.9062373528 | 2.6888 | 1.3444 |
| 640 | 4 | task | 0.9534718392 | 2.5556 | 0.6389 |
| 640 | 1 | pipeline | 0.8447707936 | 2.9070 | 2.9070 |
| 640 | 2 | pipeline | 0.8468019650 | 2.9000 | 1.4500 |
| 640 | 4 | pipeline | 0.9090972242 | 2.7013 | 0.6753 |

Ускорение рассчитано относительно SEQ baseline: `2.4367350032` для `task` и `2.4557342674` для `pipeline`. По
этим данным TBB быстрее последовательной версии, но эффективность снижается при росте числа workers из-за
накладных расходов runtime и ограниченной гранулярности блочного разбиения.

## 9. Выводы

TBB-реализация хорошо соответствует блочной природе алгоритма: основная работа описывается двумерным
диапазоном блоков результата. Главные факторы эффективности — размер `q`, стоимость одного блока и накладные
расходы runtime на разбиение диапазона.
