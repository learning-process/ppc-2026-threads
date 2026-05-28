# Сортировка Шелла с четно-нечетным слиянием Бэтчера — STL threads

- **Student:** Сахаров Александр Владимирович, группа 3823Б1ФИ3
- **Technology:** STL threads
- **Variant:** 16

## 1. Контекст

STL-версия использует стандартные потоки C++ и демонстрирует ручное управление параллелизмом. В отличие от OpenMP и
oneTBB, здесь явно создаются объекты `std::thread`, распределяются диапазоны работы и выполняется синхронизация через
`join`.

## 2. Постановка задачи

На вход поступает `std::vector<int>`, на выходе должен быть отсортированный `std::vector<int>`. Результат должен
совпадать с SEQ-версией и эталонной сортировкой.

## 3. Базовый алгоритм

Массив разбивается на чанки. Каждый чанк сортируется независимо, после чего отсортированные диапазоны объединяются
деревом попарных слияний. Слияние выполняется по уровням: сначала пары соседних чанков, затем пары уже объединенных
диапазонов и так далее.

## 4. Схема распараллеливания

На этапе сортировки для каждого чанка создается поток `std::thread`. На этапе слияния операции текущего уровня
распределяются между ограниченным числом рабочих потоков. После запуска всех потоков текущего этапа выполняется `join`,
затем буферы меняются местами.

Каждый поток пишет только в свой участок выходного буфера. Поэтому гонок данных при сортировке и слиянии не возникает, а
`mutex` и `atomic` не требуются.

## 5. Детали реализации

Файлы реализации:

- `stl/include/ops_stl.hpp`;
- `stl/src/ops_stl.cpp`.

Класс `SakharovAShellButcherSTL` возвращает тип задачи `ppc::task::TypeOfTask::kSTL`.

Минимальный листинг из `stl/src/ops_stl.cpp`:

```cpp
std::vector<std::thread> threads;
threads.reserve(worker_count);

for (std::size_t worker = 0; worker < worker_count; ++worker) {
  threads.emplace_back(MergePassThreadWorker, std::cref(source), std::ref(destination),
                       std::cref(bounds), width, current, next);
}

for (auto &thread : threads) {
  thread.join();
}
```

Фрагмент показывает ручное создание потоков и обязательное ожидание их завершения перед переходом к следующему уровню
слияния.

## 6. Проверка корректности

Корректность проверяется общим набором функциональных тестов. Эталон строится через `std::ranges::sort`. По результатам
локального запуска все 5 функциональных тестов STL прошли успешно.

## 7. Экспериментальная среда

- **CPU:** Intel Core i5-12400F, 6 cores / 12 threads, 2.50 GHz.
- **RAM:** 32 ГБ.
- **OS:** Windows 10 + WSL, Ubuntu 24.04.3 LTS.
- **IDE:** Visual Studio Code.
- **Compiler:** `g++ 13.3.0`.
- **Build system:** `CMake 3.28.3`.
- **VCS:** Git.
- **Build type:** Release-сборка проекта.
- **Perf modes:** `pipeline` и `task_run`.
- **Threads:** в таблице приведены результаты для `PPC_NUM_THREADS=2`, `4` и `8`.

## 8. Результаты

В таблице приведены средние результаты perf-тестов. Для `2` и `8` потоков среднее рассчитано по трем запускам, для `4`
потоков — по двум запускам.

| Режим | Workers | Время, с | Ускорение относительно SEQ | Эффективность |
| --- | ---: | ---: | ---: | ---: |
| pipeline | 2 | 0.0256135146 | 1.75 | 0.88 |
| task_run | 2 | 0.0233078480 | 1.94 | 0.97 |
| pipeline | 4 | 0.0140028239 | 3.20 | 0.80 |
| task_run | 4 | 0.0143011809 | 3.15 | 0.79 |
| pipeline | 8 | 0.0106910864 | 4.19 | 0.52 |
| task_run | 8 | 0.0100878557 | 4.47 | 0.56 |

## 9. Репродуцируемость

Функциональные тесты:

```bash
PPC_NUM_THREADS=8 ./build/bin/ppc_func_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_stl*'
```

Perf-тесты:

```bash
PPC_NUM_THREADS=8 ./build/bin/ppc_perf_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_stl*'
```

## 10. Выводы

STL-версия ускоряет вычисления относительно SEQ, но уступает OpenMP и oneTBB. Основная причина — ручное создание и
синхронизация потоков, которые дают больший overhead по сравнению с runtime-библиотеками.

## 11. Источники

- Лекции Сысоева А. В. по курсу «Параллельное программирование для систем с общей памятью».
- cppreference: [`std::thread`](https://en.cppreference.com/w/cpp/thread/thread), [`std::ref`,
  `std::cref`](https://en.cppreference.com/w/cpp/utility/functional/ref).
- cppreference: [`std::sort`](https://en.cppreference.com/w/cpp/algorithm/sort),
  [`std::merge`](https://en.cppreference.com/w/cpp/algorithm/merge).
