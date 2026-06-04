# Сортировка Хоара с простым слиянием

- Студент: Точилин Евгений Дмитриевич
- Группа: 3823Б1ПР3
- Вариант: 13
- Local reports: `seq/report.md`, `omp/report.md`, `tbb/report.md`, `stl/report.md`, `all/report.md`

## 1. Введение

В задаче рассматривается сортировка целочисленного массива по схеме «сортировка Хоара с простым слиянием»: вход делится
на две части, каждая часть сортируется отдельно, после чего выполняется линейное слияние отсортированных половин. Эта
задача удобна для сравнения разных моделей параллелизма, потому что в ней есть как независимые ветви quicksort, так и
неизбежный последовательный этап объединения результата.

## 2. Единая постановка задачи

Вход: `InType = std::vector<int>`, то есть одномерный массив целых чисел.

Выход: `OutType = std::vector<int>` той же длины, содержащий те же элементы в неубывающем порядке.

Ограничения:

- пустой вход считается невалидным;
- массив длины `1` допустим;
- допустимы уже отсортированные данные, обратный порядок, повторяющиеся значения и детерминированно сгенерированные
  случайные входы.

Критерий корректности:

- результат должен совпадать с эталонной сортировкой `std::ranges::sort(reference)`;
- итоговый массив должен удовлетворять `std::ranges::is_sorted(...)`.

## 3. Единая методика эксперимента

Окружение:

- CPU: AMD Ryzen 5 7500f
- Аппаратные потоки: 12
- OS: Microsoft Windows 11
- Compiler: MSVC
- Generator: Visual Studio 17 2022
- Build type: `Release`

Переменные окружения:

- `PPC_NUM_THREADS` задаёт число потоков для thread-backend-ов;
- для `ALL` дополнительно используется `PPC_NUM_PROC` как число MPI-процессов;
- для OpenMP runner синхронизирует `PPC_NUM_THREADS` и `OMP_NUM_THREADS`.

Генерация входных данных:

- functional-тесты используют детерминированные сценарии `OneElement`, `TwoElements`, `EightElements`, `RandomSize`,
  `MediumSize`, `AlreadySorted`, `ReverseSorted`, `LargeSize`, `LongSize`, `LongLongSize`;
- performance-тест использует массив размера `N = 2_000_000`, сгенерированный через `std::seed_seq{N, 0x5A17C3D2U}` и
  `std::uniform_int_distribution(-10000, 10000)`.

Сравниваемые размеры и конфигурации:

- `SEQ`: `p = 1`;
- `OMP`, `TBB`, `STL`: `1, 2, 4, 8, 12, 16` потоков;
- `ALL`: конфигурации `1x1`, `2x2`, `4x4`, `8x8`, `12x12`, `16x16`, где первая величина - число rank-ов, вторая - число
  потоков внутри rank-а.

Метрика:

- `speedup = T_seq / T_backend`;
- для `OMP`, `TBB`, `STL` эффективность считалась как `speedup / workers`;
- для `ALL` эффективность считалась как `speedup / total_workers`, где `total_workers = ranks * threads_per_rank`.

Агрегация результатов:

- для каждой конфигурации выполнялся `1` прогревочный запуск и `5` рабочих запусков;
- в сравнительные таблицы занесена медиана `task_run`.

Замечание:

- инфраструктура курса содержит и `task`, и `pipeline` режимы performance-тестов;
- в данном отчёте сводная сравнительная таблица построена по `task_run`, поскольку именно task_run замеряет
работу алгоритма RunImpl.

## 4. Сводка корректности

Все backend-ы сравнивались с одной и той же эталонной проверкой: для каждого functional-сценария строилась копия входа,
сортировалась через `std::ranges::sort(reference)`, после чего результат backend-а сравнивался с `reference`.

Использованные классы тестов:

- минимальные размеры: `OneElement`, `TwoElements`;
- малые и средние размеры: `EightElements`, `RandomSize`, `MediumSize`, `LargeSize`;
- специальные случаи: `AlreadySorted`, `ReverseSorted`;
- крупные входы: `LongSize`, `LongLongSize`.

Ограничения применимости:

- все версии отклоняют пустой вход в `ValidationImpl`;
- `ALL` требует MPI-конфигурацию и использует согласованный глобальный вход между rank-ами;
- сравнительные выводы по производительности относятся к данной машине и данной конфигурации Windows/MSVC.

## 5. Агрегированные результаты

Единый baseline:

- `SEQ`, `workers = 1`, `N = 2_000_000`, `median time = 0.0619357800 s`

Сводная таблица по `task_run`:

| backend | mode | size | workers / ranks × threads | ranks | threads_per_rank | total_workers | median time, s | speedup vs seq | efficiency | notes |
| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| SEQ | task_run | 2,000,000 | 1 | 1 | 1 | 1 | 0.0619357800 | 1.00 | 1.00 | baseline |
| OMP | task_run | 2,000,000 | 1 | - | - | 1 | 0.0645099000 | 0.96 | 0.96 | task-based quicksort |
| OMP | task_run | 2,000,000 | 2 | - | - | 2 | 0.0367070600 | 1.69 | 0.84 | устойчивый выигрыш |
| OMP | task_run | 2,000,000 | 4 | - | - | 4 | 0.0275694200 | 2.25 | 0.56 | основной выигрыш уже достигнут |
| OMP | task_run | 2,000,000 | 8 | - | - | 8 | 0.0262974800 | 2.35 | 0.29 | рост замедляется |
| OMP | task_run | 2,000,000 | 12 | - | - | 12 | 0.0256463200 | 2.41 | 0.20 | лучшая медиана OMP |
| OMP | task_run | 2,000,000 | 16 | - | - | 16 | 0.0261650000 | 2.37 | 0.15 | overhead task runtime |
| TBB | task_run | 2,000,000 | 1 | - | - | 1 | 0.0612916800 | 1.01 | 1.01 | почти baseline |
| TBB | task_run | 2,000,000 | 2 | - | - | 2 | 0.0369139000 | 1.68 | 0.84 | близко к OMP |
| TBB | task_run | 2,000,000 | 4 | - | - | 4 | 0.0279692000 | 2.21 | 0.55 | стабильное масштабирование |
| TBB | task_run | 2,000,000 | 8 | - | - | 8 | 0.0217358800 | 2.85 | 0.36 | заметный выигрыш runtime |
| TBB | task_run | 2,000,000 | 12 | - | - | 12 | 0.0207893600 | 2.98 | 0.25 | высокая эффективность по времени |
| TBB | task_run | 2,000,000 | 16 | - | - | 16 | 0.0198021800 | 3.13 | 0.20 | лучшая медиана среди thread-backend-ов |
| STL | task_run | 2,000,000 | 1 | - | - | 1 | 0.0598262600 | 1.04 | 1.04 | ручное управление потоками |
| STL | task_run | 2,000,000 | 2 | - | - | 2 | 0.0387823600 | 1.60 | 0.80 | хороший ранний выигрыш |
| STL | task_run | 2,000,000 | 4 | - | - | 4 | 0.0294734200 | 2.10 | 0.52 | нормальное масштабирование |
| STL | task_run | 2,000,000 | 8 | - | - | 8 | 0.0266623200 | 2.32 | 0.29 | join/create уже заметны |
| STL | task_run | 2,000,000 | 12 | - | - | 12 | 0.0241190000 | 2.57 | 0.21 | лучшая медиана STL |
| STL | task_run | 2,000,000 | 16 | - | - | 16 | 0.0253136400 | 2.45 | 0.15 | накладные расходы растут |
| ALL | task_run | 2,000,000 | 1 × 1 | 1 | 1 | 1 | 0.0555218800 | 1.12 | 1.12 | MPI + OMP без лишней коммуникации |
| ALL | task_run | 2,000,000 | 2 × 2 | 2 | 2 | 4 | 0.0372160800 | 1.66 | 0.42 | удачная умеренная конфигурация |
| ALL | task_run | 2,000,000 | 4 × 4 | 4 | 4 | 16 | 0.0334545600 | 1.85 | 0.12 | лучшая медиана ALL |
| ALL | task_run | 2,000,000 | 8 × 8 | 8 | 8 | 64 | 0.0985665000 | 0.63 | 0.01 | коммуникация перевешивает |
| ALL | task_run | 2,000,000 | 12 × 12 | 12 | 12 | 144 | 0.2151764200 | 0.29 | 0.00 | дорогие обмены и merge |
| ALL | task_run | 2,000,000 | 16 × 16 | 16 | 16 | 256 | 0.2831757200 | 0.22 | 0.00 | MPI-overhead доминирует |

## 6. Интерпретация различий

`SEQ` показывает стоимость алгоритма без накладных расходов runtime и служит честной базой для расчёта ускорения.

`OMP` даёт уверенный выигрыш уже на `2–4` потоках за счёт task-based распараллеливания ветвей quicksort. Дальше рост
замедляется: глубина `depth_limit = 3` ограничивает дальнейшее дробление работы, а стоимость управления задачами
перестаёт быть пренебрежимой.

`TBB` показывает лучшую медиану среди thread-backend-ов. Реализация сочетает
`parallel_invoke` на верхнем уровне, bulk-параллельную обработку фронта диапазонов через `parallel_for` и ограничение
конкуренции через `global_control`. По сути runtime oneTBB лучше балансирует множество независимых диапазонов quicksort.

`STL` оказывается близкой к `OMP`, несмотря на ручное управление потоками. Это показывает, что схема с явным `create ->
work -> join`, локальными буферами диапазонов и одним итоговым слиянием работает эффективно, но при росте числа потоков
цена ручного управления начинает сказываться.

`ALL` после улучшения дерева межпроцессных merge-ов выигрывает на `2x2` и `4x4`, но затем резко теряет эффективность.
Это не провал самой идеи гибридности, а ожидаемое следствие того, что:

- объём локальной работы на rank уменьшается;
- цена MPI-обменов и локальных merge-ов после `Recv` растёт;
- финальный `Bcast` становится заметным на больших конфигурациях.

Интересная деталь сравнения: `OMP` и `TBB` близки на малом числе workers, но затем `TBB` начинает уходить вперёд. Это
объясняется тем, что на больших входах динамическая балансировка диапазонов quicksort через oneTBB оказывается
выгоднее, чем OpenMP task runtime с фиксированной глубиной породжения задач.

## 7. Репродуцируемость

Сборка:

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
```

Functional-тесты:

```powershell
.\build\bin\ppc_func_tests.exe --gtest_filter="*tochilin_e_hoar_sort_sim_mer_seq_enabled*"

$env:PPC_NUM_THREADS='8'
.\build\bin\ppc_func_tests.exe --gtest_filter="*tochilin_e_hoar_sort_sim_mer_omp_enabled*"
.\build\bin\ppc_func_tests.exe --gtest_filter="*tochilin_e_hoar_sort_sim_mer_tbb_enabled*"
.\build\bin\ppc_func_tests.exe --gtest_filter="*tochilin_e_hoar_sort_sim_mer_stl_enabled*"

$env:PPC_NUM_THREADS='4'
$env:PPC_NUM_PROC='4'
mpiexec -n 4 .\build\bin\ppc_func_tests.exe --gtest_filter="*tochilin_e_hoar_sort_sim_mer_all_enabled*"
```

Основные performance-замеры:

```powershell
.\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_tochilin_e_hoar_sort_sim_mer_seq_enabled*"

$env:PPC_NUM_THREADS='8'
.\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_tochilin_e_hoar_sort_sim_mer_omp_enabled*"
.\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_tochilin_e_hoar_sort_sim_mer_tbb_enabled*"
.\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_tochilin_e_hoar_sort_sim_mer_stl_enabled*"

$env:PPC_NUM_THREADS='4'
$env:PPC_NUM_PROC='4'
mpiexec -n 4 .\build\bin\ppc_perf_tests.exe --gtest_filter="*task_run_tochilin_e_hoar_sort_sim_mer_all_enabled*"
```

## 8. Заключение

Для данной задачи лучшей версией по времени оказалась `TBB`, поскольку она показала наименьшую медиану и наибольший
speedup на основном размере `N = 2_000_000`. `OMP` и `STL` тоже дали устойчивое ускорение и методически
сильные реализациями той же идеи. Гибридная `ALL` после улучшения межпроцессного merge-этапа стала корректной и
полезной для умеренных конфигураций, но на больших `ranks × threads` упирается в стоимость коммуникации.

Ограничения сравнения:

- результаты относятся к одной конкретной машине и конфигурации Windows/MSVC;
- в сводной таблице отражён `task_run`;

Возможные направления развития:

- более тонкая настройка порога последовательной обработки диапазонов в `TBB` и `STL`;
- адаптивная глубина порождения задач в `OMP`;
- дальнейшее снижение коммуникационных издержек в `ALL`.

## 9. Источники

1. Документация курса <https://learning-process.github.io/parallel_programming_course/ru/index.html>.
2. OpenMP API Specification.
3. oneTBB Documentation.
4. MPI Forum Documentation.
5. Microsoft C++ Standard Library documentation for `std::thread`.
6. LLVM documentation for sanitizers and profiling tools.
7. Лекции преподавателей.

## 10. Приложение

Короткие листинги и подробности по каждой реализации вынесены в локальные отчёты:

- `seq/report.md`
- `omp/report.md`
- `tbb/report.md`
- `stl/report.md`
- `all/report.md`
