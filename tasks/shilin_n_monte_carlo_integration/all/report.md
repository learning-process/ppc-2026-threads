# Многомерное интегрирование Монте-Карло — ALL (MPI + OpenMP)

- **Student:** Шилин Никита Дмитриевич, группа 3823Б1ПР1
- **Technology:** ALL = **MPI + OpenMP**
- **Variant:** 12

---

## 1. Контекст

Версия `ALL` — **гибридная**: данные распределяются между MPI-процессами,
а внутри каждого процесса работа делится между потоками OpenMP. Такая
комбинация — самая «бюджетная» из возможных гибридов («MPI + любая
технология на выбор»; OpenMP — самый лёгкий, как и подсказывал
преподаватель), и для embarrassingly parallel Монте-Карло она даёт
наилучшее соотношение «понятность кода / выигрыш».

Внутри одного файла [`all/src/ops_all.cpp`](src/ops_all.cpp) собраны:
ветвление по рангу в `ValidationImpl`, MPI-коллектив `Allreduce`, OMP-регион
с `reduction(+:local_sum)`, барьер `MPI_Barrier`. Никаких других технологий
параллелизма в ветке `ALL` не используется — `MPI + OpenMP` хватает.

## 2. Постановка задачи

Та же, что в [`seq/report.md`](../seq/report.md). Тип задачи —
`TypeOfTask::kALL`; класс — `ShilinNMonteCarloIntegrationALL`.

В версии `ALL` валидация выполнятся **только на ранге 0**, потому что
именно туда инфраструктура курса доставляет вход; остальные ранги
пропускают проверку и сразу возвращают `true`.

```19:42:tasks/shilin_n_monte_carlo_integration/all/src/ops_all.cpp
bool ShilinNMonteCarloIntegrationALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }
  const auto &[lower, upper, n, func_type] = GetInput();
  ...
}
```

## 3. Базовый алгоритм

`RunImpl` использует **двухуровневую** декомпозицию.

- **MPI-уровень.** Глобальные индексы выборок \(0,1,\dots,n-1\) разбиваются
  по рангам **циклически**: ранг \(r\) обрабатывает все \(i\) такие, что
  \(i \equiv r \pmod{P}\), где \(P = \mathrm{num\_ranks}\). Локальное
  число итераций — \(\lceil (n - r)/P \rceil\) (формула `local_count` ниже).
  Это устраняет необходимость в `MPI_Scatterv`/`MPI_Gatherv`: исходных
  «массивов данных» не существует, выборки **детерминированно**
  генерируются по индексу через последовательность Кронекера.
- **OpenMP-уровень.** Внутри ранга цикл по локальному индексу `k`
  параллелится `#pragma omp parallel for schedule(static)` с
  `reduction(+:local_sum)`; глобальный индекс восстанавливается как
  \(i = r + k\,P\). Координата точки строится в замкнутом виде через
  `std::floor` (как в OMP-версии).
- **Глобальная редукция.** Локальные суммы складываются через
  `MPI_Allreduce(MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD)`. После этого
  **на каждом ранге** значение `global_sum` одинаково; объём
  \(V = \prod_d (\mathrm{upper}_d-\mathrm{lower}_d)\) вычисляется
  одинаково на всех рангах из локально хранящихся `lower_bounds_/upper_bounds_`.
- **Финальный барьер.** `MPI_Barrier(MPI_COMM_WORLD)` перед `return`
  фиксирует, что все ранги покинули `RunImpl` синхронно — это согласовано
  со стилем других ALL-задач в репозитории.

## 4. Межпроцессная схема (MPI)

```53:108:tasks/shilin_n_monte_carlo_integration/all/src/ops_all.cpp
bool ShilinNMonteCarloIntegrationALL::RunImpl() {
  int rank = 0;
  int num_ranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  auto dimensions = static_cast<int>(lower_bounds_.size());

  const std::vector<double> alpha = { /* sqrt(2)..sqrt(29) */ };

  int local_count = 0;
  if (rank < num_points_) {
    local_count = (((num_points_ - rank - 1) / num_ranks) + 1);
  }

  double local_sum = 0.0;

  // MSVC OpenMP does not allow non-static data members in data-sharing clauses; use `this`.
  auto *self = this;
#pragma omp parallel default(none) shared(dimensions, alpha, self, rank, num_ranks, local_count) \
    reduction(+ : local_sum)
  {
    std::vector<double> point(dimensions);
#pragma omp for schedule(static)
    for (int k = 0; k < local_count; ++k) {
      int i = rank + (k * num_ranks);
      ...
      local_sum += IntegrandFunction::Evaluate(self->func_type_, point);
    }
  }

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double volume = 1.0;
  for (int di = 0; di < dimensions; ++di) {
    volume *= (upper_bounds_[di] - lower_bounds_[di]);
  }

  GetOutput() = volume * global_sum / static_cast<double>(num_points_);
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}
```

Используются ровно две MPI-операции:

- **`MPI_Allreduce`** объединяет частные `local_sum` в `global_sum` и
  раздаёт результат всем рангам (стандартное поведение `Allreduce` по
  спецификации MPI). Это даёт согласованный `GetOutput()` на всех рангах
  без дополнительной точки рассылки.
- **`MPI_Barrier`** — точка синхронизации в конце фазы вычислений; нужен,
  чтобы `PostProcessingImpl` стартовал на всех рангах одновременно и не
  было гонок при доступе к разделяемым ресурсам инфраструктуры курса.

`Scatterv`/`Gatherv` не нужны — это и есть основная отличительная
особенность Monte Carlo от задач сортировки или контрастирования
(`PR #871`, `PR #772`), где входные данные разбираются между рангами.

## 5. Внутрипроцессная схема (OpenMP)

В OMP-регионе:

- **`default(none)`** — все переменные перечислены явно;
- **`shared(dimensions, alpha, self, rank, num_ranks, local_count)`** —
  read-only данные, поделённые между потоками. Через `self` (типа `auto*`)
  обращаемся к полям класса (`lower_bounds_`, `upper_bounds_`, `func_type_`).
  Этот трюк нужен потому, что MSVC OpenMP не разрешает использовать
  нестатические члены прямо в `shared(...)`; здесь добавлен явный
  комментарий о причине;
- **`reduction(+ : local_sum)`** — стандартная редукция;
- **`#pragma omp for schedule(static)`** — статическое равномерное
  разбиение `[0, local_count)` по потокам; для embarrassingly parallel
  это правильный выбор;
- `std::vector<double> point(dimensions)` объявлен **внутри** блока
  `parallel`, поэтому он thread-private.

Число потоков ОМП — то же, что задано через `OMP_NUM_THREADS`
(курсовый runner экспортирует `OMP_NUM_THREADS=PPC_NUM_THREADS` для
ALL-теста под `mpirun`).

## 6. Детали реализации

- Файлы: [`all/include/ops_all.hpp`](include/ops_all.hpp),
  [`all/src/ops_all.cpp`](src/ops_all.cpp).
- В `PreProcessingImpl` каждый ранг сохраняет границы и `n` локально —
  это безопасно даже без `MPI_Bcast`, потому что курсовый каркас
  передаёт `InType` в каждый ранг при создании задачи.
- `PostProcessingImpl` — `return true`; вычисленный `GetOutput()` уже
  одинаков на всех рангах (благодаря `MPI_Allreduce`).

## 7. Проверка корректности

Под `mpirun -np 2 --oversubscribe` все 8 функциональных кейсов из
[`tests/functional/main.cpp`](../tests/functional/main.cpp) проходят
(сценарий, аналогичный CI-конфигурациям курса). Без `mpirun`
функциональные ALL-кейсы помечаются как `Skipped` через стандартный
гард `func_test_util.hpp:72` — это ожидаемо.

`MPI_Allreduce` гарантирует, что `GetOutput()` на ранге 0 и
других рангах совпадает побитово (одинаковая последовательность
сложений плавающих чисел в коммутативно-ассоциативном смысле, который
поддерживает MPI для встроенных типов и операций), поэтому общий
тест корректности (`CheckTestOutputData` сверяется с аналитическим
интегралом с погрешностью \(10V/\sqrt n\)) проходит на всех рангах.

## 8. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | Apple M4 Max, 16 ядер (12P + 4E) |
| RAM | 64 GiB |
| OS | macOS 26.3.1 |
| Компилятор | Apple clang 17.0.0 |
| MPI | Open MPI 5.0.8 |
| OpenMP | libomp 21.1.8 (Homebrew) |
| Сборка | Release, `build-local` |
| Размер | \(n = 10^7\), \([0,1]^3\), `kSumSquares` |
| Повторов | 5 запусков, медиана |

Команды замера для конкретной конфигурации `P × T`:

```bash
export PPC_NUM_THREADS=$T OMP_NUM_THREADS=$T
mpirun -np $P --oversubscribe ./build-local/bin/ppc_perf_tests \
  --gtest_filter='*shilin_n_monte_carlo_integration_all*'
```

## 9. Результаты

Базовое `T_seq(task_run) = 0,030436 c`, `T_seq(pipeline) = 0,030642 c`
(см. [`seq/report.md`](../seq/report.md)).

Ниже — медиана 5 прогонов, в скобках — общее число «работников»
\(W = P \cdot T\). Эффективность считается как \(E = S / W\).

### 9.1. Режим `task_run`

| `mpirun -np` (P) ↓ \\ `PPC_NUM_THREADS` (T) → | 1 | 2 | 4 | 8 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0,030698 (S=0,99; E=99%) | 0,015917 (S=1,91; E=96%) | 0,008378 (S=3,63; E=91%) | 0,004882 (S=6,23; E=78%) |
| 2 | 0,015691 (S=1,94; E=97%) | 0,008519 (S=3,57; E=89%) | 0,004846 (S=6,28; E=79%) | 0,003909 (S=7,79; E=49%) |
| 4 | 0,008083 (S=3,77; E=94%) | 0,005516 (S=5,52; E=69%) | 0,003860 (S=7,89; E=49%) | 0,003306 (S=9,21; E=29%) |
| 8 | 0,004225 (S=7,20; E=90%) | 0,003715 (S=8,19; E=51%) | 0,003267 (S=9,32; E=29%) | **0,003140 (S=9,69; E=15%)** |

### 9.2. Режим `pipeline`

| P ↓ \\ T → | 1 | 2 | 4 | 8 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0,030637 (S=1,00; E=100%) | 0,015816 (S=1,94; E=97%) | 0,008204 (S=3,74; E=93%) | 0,006362 (S=4,82; E=60%) |
| 2 | 0,015659 (S=1,96; E=98%) | 0,008148 (S=3,76; E=94%) | 0,004668 (S=6,56; E=82%) | 0,003519 (S=8,71; E=54%) |
| 4 | 0,008020 (S=3,82; E=96%) | 0,005267 (S=5,82; E=73%) | 0,003589 (S=8,54; E=53%) | 0,003332 (S=9,20; E=29%) |
| 8 | 0,004190 (S=7,31; E=91%) | 0,003725 (S=8,23; E=51%) | 0,003331 (S=9,20; E=29%) | **0,003264 (S=9,39; E=15%)** |

**Наблюдения.**

1. **«Диагональ ядра» 12 ≤ W ≤ 16 — оптимальная зона по `S`.**
   Максимальный наблюдаемый \(S \approx 9{,}7\) (`task_run`) при
   `P=8, T=8` и \(S \approx 9{,}4\) (`pipeline`) при `P=8, T=8`. Это
   «потолок» полезной параллельности на 12P+4E ядрах M4 Max.
2. **Лучшая точка по эффективности — `P=8, T=1`** (`E=90%` в обоих
   режимах, `S=7.2`). Это «классический» MPI-режим без потоков, и
   именно он даёт самую высокую `E`.
3. **Oversubscription заметен ярко при `P=8, T≥2`:** значения `E` падают
   ниже 30%. На 16 физических ядрах 8 рангов × 2 потока = 16 потоков —
   на пределе; уже 8 × 4 = 32 потока — чистая oversubscription, и ОС
   начинает переключать контексты, теряя время.
4. **При сбалансированных конфигурациях `(P, T)` с `W ≤ 8`** эффективность
   стабильно высокая (`E ≥ 78%`): `MPI_Allreduce` + `MPI_Barrier` дают
   пренебрежимо малый overhead на данном размере задачи.
5. **Сравнение с чистым OpenMP.** Лучший OMP — `T=12, S=9.26`. Лучший
   ALL — `P=8, T=8, S=9.69`. Гибрид даёт +5% к `S` ценой втрое более
   низкой `E`. Для одного узла **OpenMP остаётся практичнее**, как
   и подсказывал преподаватель; ALL имеет смысл, когда `n` ещё больше
   или процессы реально распределены по разным узлам.

## 10. Выводы

- Гибридная связка **MPI + OpenMP** в этой задаче работает так, как и
  ожидалось от учебной ALL-версии: на 1 MPI-процессе она почти не
  отличается от чистого OpenMP, а с ростом `P` показывает
  дополнительный выигрыш до пересечения с границей физических ядер
  (12P + 4E на M4 Max).
- Лучшая точка по `S` — `P=8, T=8` (`S=9.69`); лучшая точка по `E` —
  `P=8, T=1` (`E=90%, S=7.2`). Их выводы согласованы с теорией:
  гибрид имеет смысл при достаточном объёме локальной работы и
  быстром глобальном барьере, и оба условия в Монте-Карло выполнены.
- На фоне new TBB-версии (`S=11.56`) гибридный ALL **проигрывает по `S`**
  на одноузловой машине; это связано с тем, что oneTBB лучше работает
  с гетерогенными P/E-ядрами через work-stealing, а ALL-версия
  фиксирует разбиение работы на этапе MPI и не может «украсть» итерации
  у отстающего ранга.
- Реализация **в одном файле** `all/src/ops_all.cpp` — это сознательный
  выбор: курс требует «MPI + любая технология в одной версии», а в
  `MPI + OpenMP` нет смысла дробить логику по нескольким TU.
- Версия проходит общие функциональные тесты и совместима с
  каркасом perf-тестов курса (`BaseRunPerfTests` + `MakeAllPerfTasks`),
  то есть удовлетворяет требованию преподавателя об унификации
  тестовой инфраструктуры.

## 11. Источники

- Документация курса PPC, репозиторий `ppc-2026-threads`.
- MPI Forum: спецификация MPI, описания
  [`MPI_Comm_rank`](https://www.mpi-forum.org/),
  [`MPI_Allreduce`](https://www.mpi-forum.org/),
  [`MPI_Barrier`](https://www.mpi-forum.org/).
- OpenMP API Specification — `parallel`, `for`, `reduction`,
  `default(none)`, `schedule`.
- Open MPI 5.0.8, libomp 21.1.8 — реализационные детали (Homebrew).

## 12. Чек-лист

- [x] Описаны **оба уровня** параллелизма: MPI (циклическое разбиение,
      `Allreduce`, `Barrier`) и OpenMP (`default(none)`, `shared`,
      `reduction(+:local_sum)`, `schedule(static)`).
- [x] Объяснено, почему **нет** `Scatterv/Gatherv` (нет исходных
      массивов данных, индексы детерминированные).
- [x] Указана конфигурация `ranks × threads`; нормировка
      \(E = S / (P \cdot T)\).
- [x] Корректность подтверждена под `mpirun -np 2 --oversubscribe`,
      8 функциональных кейсов проходят.
- [x] Таблицы 4×4 для `task_run` и `pipeline` (медиана 5).
- [x] Зафиксирован эффект oversubscription при `P=8`.
- [x] Этот файл — `all/report.md`, лежит ровно в `tasks/<task>/all/`.
