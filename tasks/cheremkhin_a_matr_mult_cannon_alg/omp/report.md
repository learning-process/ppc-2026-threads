# Умножение плотных матриц. Элементы типа double. Блочная схема, алгоритм Кэннона. — OMP

- Student: Черемхин Андрей Александрович
- Technology: OMP
- Variant: 1

## 1. Контекст

OpenMP-версия переносит последовательную блочную схему на многопоточное исполнение. Параллелятся три
регулярные части: копирование входных матриц в padded-буферы, вычисление независимых блоков результата и
копирование итоговой области обратно в выходной массив.

## 2. Постановка задачи

Постановка совпадает с SEQ: для двух матриц `A` и `B` размера `n x n` требуется вычислить `C = A * B`.
Последовательная версия считается эталоном корректности и источником baseline-времени.

## 3. Базовый алгоритм

Алгоритм использует виртуальную сетку блоков `q x q`, размер блока `bs = ceil(n / q)` и расширенный размер `np
= q * bs`. Каждый блок результата `C[bi, bj]` вычисляется независимо от остальных блоков, поэтому внешний цикл
по координатам блоков можно распараллелить.

## 4. Схема распараллеливания

Число потоков берётся из `ppc::util::GetNumThreads()` и задаётся вызовом `omp_set_num_threads`. В runner-е это
значение управляется переменной `PPC_NUM_THREADS`, которая также передаётся как `OMP_NUM_THREADS`.

Используются три директивы:

- `parallel for schedule(static)` для копирования строк входных матриц;
- `parallel for collapse(2) schedule(static)` для распределения пар блоков `(bi, bj)`;
- `parallel for schedule(static)` для копирования строк результата.

Для основной области `a`, `b`, `c`, `np`, `bs`, `q` являются `shared`. Индексы циклов и локальные переменные
внутри `MulAddBlock` приватны по правилам OpenMP. `reduction`, `atomic` и `critical` не требуются, потому что
каждая итерация пары `(bi, bj)` записывает только в свой блок `C[bi, bj]`.

Ключевой фрагмент:

```cpp
#pragma omp parallel for default(none) collapse(2) schedule(static) shared(a, b, c, np, bs, q, q64)
for (std::int64_t bi = 0; bi < q64; ++bi) {
  for (std::int64_t bj = 0; bj < q64; ++bj) {
    for (std::size_t step = 0; step < q; ++step) {
      const std::size_t bk = (static_cast<std::size_t>(bi) + static_cast<std::size_t>(bj) + step) % q;
      MulAddBlock(a, b, c, np, bs, static_cast<std::size_t>(bi), bk, static_cast<std::size_t>(bj));
    }
  }
}
```

`default(none)` делает область данных явной, а `schedule(static)` подходит для равномерной стоимости блоков. В
конце каждой `parallel for`-области есть неявный барьер, который безопасен: следующая фаза должна видеть
полностью подготовленные буферы.

## 5. Детали реализации

Файлы реализации: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`.

Относительно SEQ изменены только участки с независимыми циклами. Функция `MulAddBlock` сохраняет локальный
характер записи: поток работает с блоком результата, соответствующим его итерации внешнего цикла. Поэтому
гонок на `c` нет.

## 6. Проверка корректности

Корректность проверяется сравнением с oracle из функциональных тестов и с результатом SEQ на одинаковых
входах. Набор включает размеры `1`, `2`, `3`, `4`, `7`, `10`, `15`; допуск сравнения равен `1e-7`.

Дополнительно для OMP важно запускать тесты при разном числе потоков:

```bash
cmake -S . -B build -D USE_COVERAGE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_BUILD_TYPE=Releas
cmake --build build -j --parallel
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
mpirun --oversubscribe -n 1 ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_omp*'
mpirun --oversubscribe -n 2 ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_omp*'
mpirun --oversubscribe -n 4 ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_omp*'
```

## 7. Экспериментальная среда

- OS: Linux 6.6.114.1-microsoft-standard-WSL2 x86_64;
- CPU: AMD Ryzen 5 5600, 6 cores / 12 threads;
- Compiler: GCC 13.3.0 with OpenMP;
- Build type: `Release`;
- Thread control: `PPC_NUM_THREADS`, `OMP_NUM_THREADS`.

## 8. Результаты

| size | threads | mode     | median time, s | speedup vs seq | efficiency |
| ---- | ------- | -------- | -------------- | -------------- | ---------- |
| 640  | 1       | task     | 0.8571173458   | 2.8429         | 2.8429     |
| 640  | 2       | task     | 0.8821170756   | 2.7624         | 1.3812     |
| 640  | 4       | task     | 0.9496027210   | 2.5661         | 0.6415     |
| 640  | 1       | pipeline | 0.851617316    | 2.8836         | 2.8836     |
| 640  | 2       | pipeline | 0.8912803404   | 2.7553         | 1.3776     |
| 640  | 4       | pipeline | 0.9369923256   | 2.6209         | 0.6552     |

Ускорение рассчитано относительно SEQ baseline: `2.4367350032` для `task` и `2.4557342674` для `pipeline`.

## 9. Выводы

OpenMP-реализация является прямым и компактным способом распараллелить блочный алгоритм. Она сохраняет
структуру SEQ, не требует ручной синхронизации и должна быть эффективной на достаточно больших матрицах, где
стоимость вычисления блоков превышает накладные расходы создания parallel regions.
