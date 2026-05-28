# Сортировка Хоара с нечётно-чётным слиянием Бэтчера - OMP

- **Студент:** Сметанин Дмитрий Владимирович, группа 3823Б1ПР3
- **Технология:** OMP (OpenMP)
- **Вариант:** 14

---

## 1. Контекст

Параллельная версия той же сортировки, что в SEQ: при длинах подотрезка выше порога используется
рекурсивное разбиение с **задачами OpenMP** (`omp task` / `taskwait`); для коротких отрезков
вызывается тот же итеративный `HoarSortBatcherSeq`, что и в последовательной реализации.

---

## 2. Постановка задачи

**Вход / выход.** `std::vector<int>` - отсортированный по неубыванию `std::vector<int>` той же длины.

Краткая постановка и описание odd-even merge после разбиения Хоара - в SEQ.

---

## 3. Базовый алгоритм

Последовательное ядро совпадает с SEQ: разбиение Хоара, затем `OddEvenMerge` на отрезке, стековая
обработка в `HoarSortBatcherSeq` для малых кусков.

---

## 4. Схема распараллеливания

**Порог.** `kTaskCutoff = 1000`: если `hi - lo < kTaskCutoff`, выполняется только `HoarSortBatcherSeq`
без порождения задач.

**Область параллелизма.** Один регион `#pragma omp parallel default(none) shared(data, n)` с
`#pragma omp single` на корне (`RunImpl` передаёт ссылку на выходной вектор как `data` и его размер
`n`); рекурсивная функция `HoarSortBatcherOMPImpl` порождает две подзадачи:

```cpp
// File: omp/src/ops_omp.cpp
#pragma omp task default(none) shared(arr) firstprivate(lo, p)
HoarSortBatcherOMPImpl(arr, lo, p);
#pragma omp task default(none) shared(arr) firstprivate(hi, p)
HoarSortBatcherOMPImpl(arr, p + 1, hi);
#pragma omp taskwait
```

**Атрибуты данных.** Массив `arr` - **shared**; границы `lo`, `hi`, `p` - **firstprivate**, чтобы каждая
задача получила свои копии границ.

**Барьеры.** После пары `omp task` следует **`taskwait`** - дочерние задачи должны завершиться до
возврата из текущего вызова.

**Число потоков.** Задаётся через `PPC_NUM_THREADS` (runner дублирует в `OMP_NUM_THREADS`); в коде нет
`omp_set_num_threads`.

---

## 5. Детали реализации

**Файлы.** `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp` - класс `SmetaninDHoarSortOMP`.

**Отличия от SEQ.** В `RunImpl` при `n > 1` вызывается `HoarSortBatcherOMPImpl` внутри `parallel` +
`single`; иначе поведение пайплайна как у SEQ.

**Гонки.** Запись в разные подотрезоны после разбиения Хоара не пересекается между задачами;
`OddEvenMerge` вызывается в каждой задаче на своём текущем отрезке до порождения потомков — порядок
согласован с логикой алгоритма.

---

## 6. Проверка корректности

Общие тесты: `tests/functional/main.cpp` (`SmetaninDRunFuncTests`),
`tests/performance/main.cpp` (`SmetaninDRunPerfTests`). Результат сравнивается с эталоном сортировки
(функционально) или проверяется `is_sorted` (перф).

**Сравнение с SEQ.** При одинаковом входе выход должен совпадать с SEQ для всех размеров из
`kTestParam` и для перф-входа длины **1 000 000**.

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

**Baseline.** Время SEQ на перф-тесте: **0,493 с** (см. в SEQ); режим `task_run`.

**Что такое «эффективность».** `efficiency = (T_seq / T_omp) / threads · 100%`.

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

- **1 поток:** 0,512 с; ускорение 0,96; эффективность 96,2%.
- **2 потока:** 0,299 с; ускорение 1,65; эффективность 82,4%.
- **4 потока:** 0,258 с; ускорение 1,91; эффективность 47,7%.
- **8 потоков:** 1,023 с; ускорение 0,48; эффективность 6,0%.

*(Режим `task_run`, среднее по 5 прогонам; ускорение и эффективность относительно **T_seq = 0,493 с**.)*

При **8 потоках** на тестовой машине зафиксировано **замедление** относительно эталона — типичный
признак перегруза планировщика и конкуренции за кэш при избыточном параллелизме для данной задачи.

---

## 9. Выводы

OpenMP-задачи дают ускорение около **×1,65–1,91** при 2 и 4 потоках; при одном потоке время близко к
SEQ с небольшим штрафом на обслуживание региона `parallel`/`single`. На **8 потоках** измерение
показывает деградацию — имеет смысл ограничивать `OMP_NUM_THREADS` или повышать порог `kTaskCutoff`
под железо.
