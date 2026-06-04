# Сортировка Хоара с нечётно-чётным слиянием Бэтчера - TBB

- **Студент:** Сметанин Дмитрий Владимирович, группа 3823Б1ПР3
- **Технология:** TBB (oneTBB)
- **Вариант:** 14

---

## 1. Контекст

Версия использует **волновую** обработку набора отрезков: на каждой итерации для всех текущих
интервалов `[l, r]` выполняется `tbb::parallel_for` по индексам отрезков; большие отрезки делятся
разбиением Хоара и odd-even merge и попадают в следующую «волну» через `tbb::concurrent_vector`.
Малые отрезки обрабатываются последовательным `HoarSortBatcherSeq` (как в SEQ).

---

## 2. Постановка задачи

Вход / выход - `std::vector<int>`, результат - неубывающая сортировка. Подробности постановки - в SEQ.

---

## 3. Базовый алгоритм

Тот же последовательный кусок (`HoarePartition`, `OddEvenMerge`, стековый `HoarSortBatcherSeq`), что и
в SEQ, используется как база для листьев и малых диапазонов.

---

## 4. Схема распараллеливания

**Примитив.** `oneapi/tbb::parallel_for` по диапазону индексов `[0, cur.size())`, где `cur` -
`std::vector<std::pair<int,int>>` активных отрезков на текущей волне.

**Сбор следующей волны.** `tbb::concurrent_vector<std::pair<int,int>> next`: каждая задача может
добавить до двух новых пар границ `{l, p}` и `{p+1, r}` после разбиения и merge.

**Порог.** `kTaskCutoff = 1000`: при `r - l < kTaskCutoff` вызывается только `HoarSortBatcherSeq(arr,
l, r)` без расщепления волны.

**Разбиение работы.** Индексное `parallel_for` равномерно распределяет **разные отрезки** между
потоками; размер чанка для планировщика задаёт рантайм TBB по умолчанию (явный `blocked_range` с
grainsize в коде не задан).

**Конкуренция.** Ограничение числа потоков задаётся окружением курса (`PPC_NUM_THREADS`), без явного
`global_control` в реализации.

**Фрагмент:**

```cpp
// File: tbb/src/ops_tbb.cpp
tbb::parallel_for(static_cast<std::size_t>(0), cur_sz, [&](std::size_t idx) {
  const int l = cur[idx].first;
  const int r = cur[idx].second;
  // порог -> HoarSortBatcherSeq; иначе partition, OddEvenMerge, next.push_back(...)
});
cur.assign(next.begin(), next.end());
```

---

## 5. Детали реализации

**Файлы.** `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp` - класс `SmetaninDHoarSortTBB`.

**Зависимости.** Подключены `oneapi/tbb/parallel_for.h` и `oneapi/tbb/concurrent_vector.h`.

**Корректность памяти.** Разные задачи на одной волне обрабатывают непересекающиеся отрезки до
порождения потомков; границы следующей волны формируются после локального разбиения на текущем отрезке.

---

## 6. Проверка корректности

Те же тесты, что для SEQ/OMP: `SmetaninDRunFuncTests`, `SmetaninDRunPerfTests`. Сравнение с SEQ по
выходу на всех функциональных размерах; на перф-тесте — `is_sorted`.

---

## 7. Экспериментальная среда

- **Процессор:** AMD Ryzen 5 5600H with Radeon Graphics
- **ОЗУ:** 16 ГБ
- **ОС:** Microsoft Windows 10
- **Компилятор:** MSVC 14.43 (Visual Studio 2022 Community)
- **Сборка:** Release
- **Open MPI:** MS-MPI (`mpiexec`)
- **oneTBB:** 2022.3.0
- **OpenMP:** OpenMP в составе MSVC

**Baseline.** Время SEQ на перф-тесте: **0,493 с** (`task_run`; см. в SEQ).

**Эффективность.** `efficiency = (T_seq / T_tbb) / threads · 100%`.

**Команды запуска.**

```bash
git submodule update --init --recursive --depth=1
cmake -S . -B build -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

export PPC_NUM_THREADS=4
scripts/run_tests.py --running-type=threads --counts 1 2 4
scripts/run_tests.py --running-type=performance
```

---

## 8. Результаты

- **1 поток:** 0,503 с; ускорение 0,98; эффективность 97,9%.
- **2 потока:** 0,298 с; ускорение 1,66; эффективность 82,8%.
- **4 потока:** 0,299 с; ускорение 1,65; эффективность 41,2%.
- **8 потоков:** 1,114 с; ускорение 0,44; эффективность 5,5%.

*(Режим `task_run`, среднее по 5 прогонам; ускорение от **T_seq = 0,493 с**.)*

При 2 потоках результат сопоставим с OMP; при 4 прирост почти не растёт за счёт структуры волн и
накладных расходов. При **8 потоках** — сильное замедление (аналогично OMP).

---

## 9. Выводы

TBB даёт около **×1,65** ускорения при 2–4 потоках на этом железе; рост до 8 потоков на ноутбуке
оказался контрпродуктивным. Волновая схема с `concurrent_vector` добавляет синхронизацию между
уровнями разбиения.
