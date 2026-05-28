# Вычисление многомерных интегралов методом прямоугольников — STL

- Студент: Галкин Данила Алексеевич
- Технология: STL
- Вариант: 9

## 1. Контекст

STL-реализация использует `std::thread` для ручного параллелизма без внешних библиотек.
Число потоков задаётся через `PPC_NUM_THREADS`. Разбиение диапазона ячеек — статическое блочное,
аналогичное MPI-разбиению в ALL-реализации.

## 2. Постановка задачи

Идентична SEQ: входные данные — функция `f`, вектор границ и `n > 0`; выход — значение
многомерного интеграла методом средней точки.

## 3. Описание алгоритма и схемы распараллеливания

### Блочное разбиение

Общий диапазон `[0, total_cells)` делится на `num_threads` блоков:

```text
chunk     = total_cells / num_threads
remainder = total_cells % num_threads

begin[tid] = tid * chunk + min(tid, remainder)
end[tid]   = begin[tid] + chunk + (tid < remainder ? 1 : 0)
```

Первые `remainder` потоков получают на одну ячейку больше — дисбаланс не превышает 1 ячейку.

### Запуск и синхронизация потоков

```cpp
// File: stl/src/ops_stl.cpp
std::vector<std::thread> threads(num_threads - 1);
for (int ti = 1; ti < num_threads; ++ti) {
    threads[ti - 1] = std::thread(worker, ti);
}
worker(0);  // главный поток выполняет долю tid=0

double sum = partial_sums[0];
for (int ti = 1; ti < num_threads; ++ti) {
    threads[ti - 1].join();
    sum += partial_sums[ti];
}
GetOutput() = sum * cell_v;
```

Все `num_threads − 1` рабочих потоков запускаются **перед** тем, как главный поток вызывает
`worker(0)`. Только после этого начинается `join()`. Такая схема гарантирует настоящий
параллелизм: все потоки работают одновременно, а не по очереди.

### Тело worker-функции

```cpp
auto worker = [&](int tid) {
    const std::int64_t begin = ...;
    const std::int64_t end   = ...;
    std::vector<double> x(dim);  // одна аллокация на поток
    double local_sum = 0.0;
    for (std::int64_t linear_idx = begin; linear_idx < end; ++linear_idx) {
        auto tmp = static_cast<std::size_t>(linear_idx);
        for (std::size_t i = 0; i < dim; ++i) {
            const std::size_t idx_i = tmp % static_cast<std::size_t>(n);
            tmp /= static_cast<std::size_t>(n);
            x[i] = borders[i].first + ((static_cast<double>(idx_i) + 0.5) * h[i]);
        }
        local_sum += func(x);
    }
    partial_sums[tid] = local_sum;  // запись в независимый слот, без mutex
};
```

Каждый поток пишет только в `partial_sums[tid]` — слоты независимы, гонок нет. Вектор `x`
создаётся один раз на поток (не на итерацию), что устраняет главный overhead OMP-реализации.

## 4. Детали реализации

**Файлы:** `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`

Вспомогательные функции (безымянное пространство имён):

- `ComputeGridParams` — вычисляет шаги `h[i]` и объём ячейки;
- `ComputeTotalCells` — проверяет переполнение и возвращает `total_cells` как `int64_t`.

Метки:

- `ValidationImpl` — проверяет корректность входа;
- `PreProcessingImpl` — сбрасывает выход в 0;
- `RunImpl` — разбивает диапазон, создаёт потоки, собирает результаты;
- `PostProcessingImpl` — проверяет `isfinite(result_)`.

## 5. Проверка корректности

Набор функциональных тестов идентичен SEQ (23 параметрических теста). Результат сравнивается
с допуском 1e-4. Все тесты проходят.

## 6. Экспериментальная среда

| Параметр | Значение |
| -------- | -------- |
| CPU | Apple M4 Pro |
| Физических ядер | 12 |
| RAM | 24 GB |
| ОС | macOS 26.3.1 |
| Компилятор | Apple clang 21.0.0 |
| Тип сборки | Release |
| Функция теста | sin(x₀)·sin(x₁) на [0,π]², n=10 000 |

Команды запуска:

```bash
cd ppc-2026-threads/local_checks/galkin_d_multidim_integrals_rectangles/build
cmake .. && make -j$(nproc)
./bin/galkin_func_tests
PPC_NUM_THREADS=4 ./bin/galkin_perf_tests --gtest_filter="*stl*"
```

## 7. Результаты

SEQ-эталон: task\_run = 0.516 с.

| Потоков | task\_run, с | Ускорение vs SEQ | Ускорение vs STL T=1 |
| ------- | ------------ | ---------------- | -------------------- |
| 1 | 1.027 | 0.50× | 1.00× |
| 2 | 0.902 | 0.57× | 1.14× |
| 4 | 0.625 | 0.83× | 1.64× |
| 8 | 0.488 | 1.06× | 2.10× |

STL при T=1 работает медленнее SEQ (~2×): накладные расходы на создание потока и блочное
разбиение добавляют overhead, которого нет в однопоточном SEQ.
При T=8 STL впервые превосходит SEQ (1.06×) благодаря распределению работы по 8 ядрам
при одной аллокации `x` на поток.

## 8. Выводы

STL-реализация правильно реализует ручной параллелизм: все потоки запускаются до вызова
`worker(0)` и `join()`, что обеспечивает истинную параллельность. Отсутствие per-iteration
аллокаций (в отличие от OMP) даёт лучшую масштабируемость. При T=8 удаётся превысить SEQ
на 6%. Для достижения сопоставимых с TBB результатов потребовалось бы задействовать все 12 ядер
и устранить overhead запуска потоков (пул потоков вместо создания «на лету»).
