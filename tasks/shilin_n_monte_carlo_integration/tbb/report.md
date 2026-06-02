# Многомерное интегрирование Монте-Карло — oneTBB

- **Student:** Шилин Никита Дмитриевич, группа 3823Б1ПР1
- **Technology:** oneTBB
- **Variant:** 12

---

## 1. Контекст

Цель — переписать главный цикл `RunImpl` через **задачно-ориентированный
примитив** oneTBB. Для embarrassingly parallel редукции естественный выбор —
`oneapi::tbb::parallel_reduce` с диапазоном `tbb::blocked_range<int>`. Этот
примитив сам строит дерево задач, разбивает диапазон, ведёт собственные
worker-потоки, поэтому в коде нет ни ручного `join`, ни `reduction`-директив —
вся семантика накопления выражена через объект-аккумулятор.

Для согласования с курсовым каркасом и для построения честной таблицы
масштабирования по `PPC_NUM_THREADS` в реализацию добавлены два
обязательных элемента (см. §4): **`tbb::global_control`** для ограничения
числа worker-ов и **явный `grainsize`** в `blocked_range`, чтобы
тяжеловесный планировщик oneTBB не дробил диапазон до итераций по 1
точке (когда стоимость одной точки — несколько FLOP).

## 2. Постановка задачи

Та же, что в [`seq/report.md`](../seq/report.md):
\(I \approx V \cdot \tfrac{1}{n}\sum_i f(x^{(i)})\). Тип задачи — `kTBB`,
класс — `ShilinNMonteCarloIntegrationTBB`.

## 3. Базовый алгоритм

Тело редукции переиспользует ту же замкнутую формулу, что и OMP
(см. [`omp/report.md`](../omp/report.md), §3): `i+1`-я точка строится через
`std::floor` без межитерационной зависимости. Это позволяет TBB свободно
делить диапазон `[0, n)` на произвольные поддиапазоны.

## 4. Схема распараллеливания

```52:89:tasks/shilin_n_monte_carlo_integration/tbb/src/ops_tbb.cpp
bool ShilinNMonteCarloIntegrationTBB::RunImpl() {
  auto dimensions = static_cast<int>(lower_bounds_.size());

  const std::vector<double> alpha = { /* sqrt(2)..sqrt(29) */ };

  // Honor PPC_NUM_THREADS so reports can construct a real T-scaling table.
  const int num_threads = std::max(1, ppc::util::GetNumThreads());
  tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism,
                                   static_cast<std::size_t>(num_threads));
  // Explicit grain prevents oneTBB from over-splitting a very light per-iteration body.
  const int grain = std::max(1024, num_points_ / (num_threads * 4));

  double sum = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, num_points_, grain), 0.0,
      [&](const tbb::blocked_range<int> &range, double local_sum) {
        std::vector<double> point(dimensions);
        for (int i = range.begin(); i < range.end(); ++i) {
          for (int di = 0; di < dimensions; ++di) {
            double val = 0.5 + (static_cast<double>(i + 1) * alpha[di]);
            double current = val - std::floor(val);
            point[di] = lower_bounds_[di] + ((upper_bounds_[di] - lower_bounds_[di]) * current);
          }
          local_sum += IntegrandFunction::Evaluate(func_type_, point);
        }
        return local_sum;
      },
      std::plus<>());
```

Расшифровка по компонентам:

- **`tbb::global_control(max_allowed_parallelism, num_threads)`** —
  ограничение числа worker-ов TBB значением `PPC_NUM_THREADS`. Без него
  oneTBB запустила бы свой default-pool (на M4 Max — 16 потоков
  независимо от настройки), и таблицы масштабирования по `T` стали бы
  бессмысленными. Установленный паттерн в курсе — см.
  `shkrebko_m_calc_of_integral_rect/tbb/src/ops_tbb.cpp:62-63`.
- **`tbb::blocked_range<int>(0, num_points_, grain)`** — итерируемый
  диапазон с **явным `grainsize`**, рассчитанным как
  `max(1024, n/(T*4))`. На `n=10^7, T=8` получается grain ≈ 312k — это
  даёт ~32 крупных задачи, что хорошо подходит и для work-stealing TBB,
  и для embarrassingly parallel тела с дешёвой итерацией. Без явного
  grain oneTBB могла бы дробить диапазон вплоть до одиночных индексов,
  что съедало весь параллелизм накладными расходами планировщика
  (это и наблюдалось в первой версии задачи, где `S ≈ 1.04`).
- **Лямбда тела** принимает `(range, local_sum)` и возвращает
  обновлённый аккумулятор. Это **функциональная** редукция, в отличие от
  OMP-`reduction`: TBB сама создаёт нужное число локальных копий
  `local_sum`, инициализирует их единичным значением `0.0`, наполняет в
  параллельных задачах и затем сворачивает через **identity-комбинатор**
  `std::plus<>()` в одно глобальное значение.
- **`std::vector<double> point(dimensions)`** объявлен **внутри** лямбды —
  свой буфер у каждой задачи, никаких общих записей и false sharing нет.

## 5. Детали реализации

- Файлы: [`tbb/include/ops_tbb.hpp`](include/ops_tbb.hpp),
  [`tbb/src/ops_tbb.cpp`](src/ops_tbb.cpp).
- Заголовки `oneapi/tbb/blocked_range.h`,
  `oneapi/tbb/global_control.h` и `oneapi/tbb/parallel_reduce.h`
  подключаются явно (стиль oneTBB 2022.x). Дополнительно подключается
  `util/include/util.hpp` ради `ppc::util::GetNumThreads()`. Собрано
  на oneTBB 2022.3.0 (Homebrew).
- `ValidationImpl` совпадает с SEQ (лежит в одном процессе, MPI не нужен).
- `PreProcessingImpl` повторяет SEQ.

## 6. Проверка корректности

8 параметризованных кейсов из
[`tests/functional/main.cpp`](../tests/functional/main.cpp) — все проходят.
Точность не зависит от числа worker-ов TBB, потому что лямбда вычисляет
`Evaluate(func_type, point)` независимо для каждого индекса \(i\).

## 7. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | Apple M4 Max, 16 ядер (12P + 4E) |
| RAM | 64 GiB |
| OS | macOS 26.3.1 |
| Компилятор | Apple clang 17.0.0 |
| oneTBB | 2022.3.0 (Homebrew) |
| Сборка | Release, `build-local` |
| Размер | \(n = 10^7\), \([0,1]^3\), `kSumSquares` |
| Повторов | 5 запусков, медиана |

Команда замера для конкретного `T`:

```bash
export PPC_NUM_THREADS=$T OMP_NUM_THREADS=$T
./build-local/bin/ppc_perf_tests --gtest_filter='*shilin_n_monte_carlo_integration_tbb*'
```

## 8. Результаты

Базовое `T_seq(task_run) = 0,030436 c`,
`T_seq(pipeline) = 0,030642 c` (см. [`seq/report.md`](../seq/report.md)).

### 8.1. Режим `task_run`

| `PPC_NUM_THREADS` | время, с | \(S\) | \(E\) |
| ---: | ---: | ---: | ---: |
| 1 | 0,029259 | 1,04 | 104,0% |
| 2 | 0,014979 | 2,03 | 101,6% |
| 4 | 0,007867 | 3,87 | 96,7% |
| 8 | 0,004286 | 7,10 | 88,7% |
| 16 | **0,002632** | **11,56** | **72,3%** |

### 8.2. Режим `pipeline`

| `PPC_NUM_THREADS` | время, с | \(S\) | \(E\) |
| ---: | ---: | ---: | ---: |
| 1 | 0,029060 | 1,05 | 105,4% |
| 2 | 0,015024 | 2,04 | 102,0% |
| 4 | 0,007819 | 3,92 | 98,0% |
| 8 | 0,004295 | 7,14 | 89,2% |
| 16 | **0,002736** | **11,20** | **70,0%** |

**Наблюдения.**

1. **TBB после правки — лучший backend по абсолютному `S`.** При `T=16`
   \(S \approx 11{,}6\) — это превосходит как OMP (`S=9.3` при `T=12`),
   так и STL (`S=4.97` при `T=16`). Причина: oneTBB с явным `grainsize`
   и `global_control` строит ровно столько крупных задач, сколько нужно
   для work-stealing на 16 ядрах, и активно перераспределяет работу
   между ядрами в течение прогона. Для гетерогенного M4 Max это
   принципиальный плюс перед OMP: «отстающие» E-ядра не блокируют всех,
   потому что P-ядра «крадут» у них работу через TBB-арена.
2. **Эффективность 72% при `T=16`** объясняется тем, что 4 E-ядра в
   ~2 раза медленнее P-ядер, и часть worker-ов в среднем простаивает
   меньше, но всё равно медленнее, чем теоретическое 100%.
3. **До `T=8` рост близок к линейному** (`E ≥ 89%`); это типичная
   картина для well-tuned TBB на embarrassingly parallel задаче.
4. `task_run` и `pipeline` дают одинаковые цифры в пределах 3% — для
   лёгкой задачи разница в каркасе незаметна.

### 8.3. Историческая справка о «дефолтном» TBB

В **первой** версии задачи в коде **не было** ни `tbb::global_control`,
ни явного `grainsize`. Это давало:

- `S ≈ 1,04` независимо от `PPC_NUM_THREADS` (oneTBB всегда брала 16
  потоков по умолчанию);
- неблагоприятное соотношение «overhead планировщика / полезная работа»,
  потому что auto-partitioner дробил диапазон до очень мелких задач.

Текущая версия (см. §4) исправляет оба недостатка — отсюда и реальные
цифры в §8.1–§8.2. Изменение зафиксировано в отдельном `[FIXED]` PR.

## 9. Выводы

- TBB-реализация **корректна** (проходит 8 функциональных тестов) и
  **самая быстрая** из четырёх потоковых backend-ов: `S = 11,56` при
  `T = 16`, эффективность 72%.
- Залог такого результата — сочетание `tbb::global_control` (TBB
  уважает `PPC_NUM_THREADS`) и явного `grainsize` в `blocked_range`
  (TBB не дробит лёгкие итерации до single-task). Без этих двух
  настроек TBB на этой задаче выглядела бы как «SEQ с обвязкой»
  (`S ≈ 1`) — что и было до правки.
- Главное преимущество TBB перед OMP на M4 Max — **work stealing**: на
  гетерогенных P/E-ядрах работа перераспределяется в рантайме, а
  «отстающие» E-ядра не блокируют барьер.

## 10. Источники

- Документация курса PPC, репозиторий `ppc-2026-threads`.
- oneAPI Threading Building Blocks (UXL Foundation), 2022.x —
  описания `parallel_reduce`, `blocked_range`, auto-partitioner,
  `global_control`.
- Homebrew: формула `tbb` версии 2022.3.0.

## 11. Чек-лист

- [x] Описаны `parallel_reduce`, `blocked_range`, identity, комбинатор.
- [x] Зафиксирован переход от «дефолтного» TBB к версии с
      `tbb::global_control` и явным `grainsize` (`[FIXED]` PR), с
      описанием эффекта (`S` с 1,04 до 11,56).
- [x] Корректность подтверждена 8 функциональными кейсами при
      `PPC_NUM_THREADS ∈ {1, 8, 16}`.
- [x] Таблицы по \(T \in \{1, 2, 4, 8, 16\}\) для `task_run` и
      `pipeline` (медиана 5).
- [x] В выводах объяснено, почему TBB здесь обходит OMP (work stealing
      на гетерогенных ядрах M4 Max).
- [x] Этот файл — `tbb/report.md`, лежит ровно в `tasks/<task>/tbb/`.
