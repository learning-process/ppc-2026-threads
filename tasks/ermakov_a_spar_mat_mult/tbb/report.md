# Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – строковый (CRS) - TBB

- Студент: Ермаков Алексей Викторович
- Группа: 3823Б1ПР3
- Технология: TBB
- Вариант: 6

## 1. Контекст

TBB-версия интересна тем, что описывает распараллеливание через задачно-
ориентированную декомпозицию диапазонов строк, а не через явное управление
потоками.

## 2. Постановка задачи

На вход подаются две комплексные CRS-матрицы `A` и `B`, на выходе требуется
получить `C = A * B` в том же формате.

Baseline-реализация из `seq/report.md` остается эталоном:

- по корректности результата;
- по входным ограничениям;
- по формуле `speedup = T_seq / T_tbb`.

## 3. Базовый алгоритм

Математическая логика та же, что и в SEQ:

1. для строки `i` матрицы `A` перебираются все `A[i, j]`;
2. по каждому такому элементу просматривается строка `j` матрицы `B`;
3. произведения `A[i, j] * B[j, k]` накапливаются по столбцам `k`;
4. затронутые столбцы сортируются;
5. строка переносится в итоговый CRS.

TBB меняет только организацию исполнения: строки делятся на диапазоны
`blocked_range`, а рабочие буферы хранятся в thread-local workspace.

## 4. Схема распараллеливания

Выбранный примитив - `oneapi::tbb::parallel_for`.

```cpp
tbb::parallel_for(tbb::blocked_range<int>(0, m, grain_size), [&](const tbb::blocked_range<int> &range) {
  auto &local = workspace.local();

  for (int i = range.begin(); i != range.end(); ++i) {
    AccumulateRowProducts(i, local.row_vals, local.row_mark, local.used_cols);
    SortUsedCols(local.used_cols);
    CollectRowValues(local.row_vals, local.used_cols, row_cols[static_cast<std::size_t>(i)],
                     row_values[static_cast<std::size_t>(i)]);
  }
});
```

Расшифровка решений:

- диапазон работы задается как `blocked_range<int>(0, m, grain_size)`;
- каждая задача получает не одну строку, а блок строк, что уменьшает overhead;
- `grain_size` выбирается функцией `ResolveGrainSize(m)` и примерно равен `m / 16`,
  но не меньше `1`;
- `partitioner` явно не задается, используется поведение oneTBB по умолчанию;
- `global_control::max_allowed_parallelism` в этой реализации не используется;
- фактическая конкуренция ограничивалась инфраструктурой запуска через
  `PPC_NUM_THREADS`.

## 5. Детали реализации

Файлы:

- `tbb/include/ops_tbb.hpp`
- `tbb/src/ops_tbb.cpp`

Ключевые элементы реализации:

- `RowWorkspace` хранит `row_vals`, `row_mark`, `used_cols`;
- `tbb::enumerable_thread_specific<RowWorkspace>` дает каждому worker-у
  собственный workspace;
- `AccumulateRowProducts` выполняет арифметическое накопление;
- `CollectRowValues` переносит строку в построчные контейнеры;
- итоговый CRS собирается последовательно после `parallel_for`.

Ложное разделение кэша через специальные allocator-ы здесь отдельно не
оптимизировалось; ставка сделана на простую и безопасную thread-local схему.

## 6. Проверка корректности

Корректность TBB-версии проверялась тем же набором тестов, что и остальные
backend-ы:

- `SmallFixed (3x3)`;
- детерминированные разреженные входы размеров `10`, `20`, `30`;
- сравнение с плотным эталоном `DenseMul`;
- сравнение с baseline-версией.

Это покрывает и граничные размеры, и разные степени разреженности.

## 7. Экспериментальная среда

- CPU: AMD Ryzen 5 1600
- RAM: 8 GB
- OS: Windows 11
- compiler: MSVC
- build type: `Release`

Способ ограничения конкуренции:

- `global_control` в коде не используется;
- фактическое число workers задавалось через `PPC_NUM_THREADS`.

Команды запуска:

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake -S . -B build -D ENABLE_ADDRESS_SANITIZER=ON -D CMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --config Release --parallel
cd build\bin
$env:PPC_NUM_THREADS='8'; .\ppc_func_tests.exe --gtest_filter="*ermakov_a_spar_mat_mult_tbb_enabled*"
$env:PPC_NUM_THREADS='8'; .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_tbb_enabled"
```

Размер performance-задачи:

- `15000 x 15000`;
- плотность `0.001`;
- детерминированный CRS-вход.

## 8. Результаты

Последовательный baseline:

- `T_seq = 0.1341834602 c`.

Таблица TBB-замеров:

| Workers | Time, s      | Speedup | Efficiency |
| ------: | -----------: | ------: | ---------: |
| 1       | 0.1847540202 | 0.726   | 0.726      |
| 2       | 0.1124937002 | 1.193   | 0.596      |
| 4       | 0.0786577202 | 1.706   | 0.426      |
| 8       | 0.0699605602 | 1.918   | 0.240      |
| 12      | 0.0697641002 | 1.923   | 0.160      |
| 16      | 0.0698530402 | 1.921   | 0.120      |

Комментарий по balancing и overhead:

- TBB показывает лучший результат среди всех shared-memory реализаций;
- диапазон `4–8` workers дает основной рост производительности;
- после `8–12` workers ускорение практически насыщается;
- выбранный `grain_size` оказался достаточно крупным, чтобы не утонуть в overhead,
  и достаточно мелким, чтобы runtime мог балансировать работу.

## 9. Выводы

TBB оказалась лучшей технологией на этой задаче. Сильная сторона реализации -
сочетание `parallel_for`, `blocked_range` и thread-local workspace. Overhead
начинает мешать тогда, когда дальнейший рост workers уже упирается не в runtime,
а в память и в структуру самого алгоритма.
