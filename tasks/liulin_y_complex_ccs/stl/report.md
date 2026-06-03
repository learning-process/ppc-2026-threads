# Умножение разреженных матриц с комплексными элементами - STL

- **Student:** Люлин Ярослав Сергеевич
- **Technology:** STL
- **Variant:** 7

## 1. Контекст

Реализация на базе стандартной библиотеки C++ (STL) подразумевает низкоуровневое
управление потоками через `std::thread`. В отличие от OMP или TBB, здесь разработчик
самостоятельно распределяет итерационное пространство и управляет жизненным циклом
вычислительных единиц. Цель работы — реализовать эффективное параллельное умножение
без использования сторонних фреймворков и библиотек спецификаций.

## 2. Постановка задачи

- **Вход/Выход:** Две матрицы CCS `A` и `B`, результат `C = A * B` в формате CCS.
- **Особенности реализации:** Использование `std::map` для автоматической сортировки
  элементов и `std::thread` для параллельных вычислений.
- **Задача:** Разделить общий объем работы (количество ненулевых элементов матрицы `A`)
  на равные части между потоками, минимизировав при этом время на синхронизацию.

## 3. Схема распараллеливания

В данной реализации применена стратегия **декомпозиции по данным (Data Decomposition)**
на уровне отдельных элементов матрицы `A`.

1. **Подготовка:** Матрица `B` транспонируется для обеспечения быстрого доступа к строкам.
   Строится вспомогательный массив индексов столбцов для элементов матрицы `A`.
2. **Разбиение нагрузки:** Общее количество ненулевых элементов матрицы `A` ($nnz\_A$)
   делится на число доступных аппаратных ядер процессора (`hardware_concurrency`).
3. **Параллельная обработка:** Каждый поток получает свой диапазон (chunk) индексов
   элементов матрицы `A`.
4. **Изоляция данных:** Для исключения состояния гонки (Data Race) и отказа от блокировок
   (`mutex`), каждый поток владеет собственной локальной структурой данных
   для накопления сумм — `LocalDict` (на базе `std::map`).
5. **Сборка:** После завершения работы всех потоков (`join()`), мастер-поток объединяет
   локальные словари в один глобальный и формирует финальную структуру CCS.

## 4. Особенности реализации

### Использование `std::map` как аккумулятора

Ключом в словаре является пара `std::pair<int, int>`, представляющая координаты
`{столбец, строка}`. Выбор именно такого порядка в паре обусловлен тем, что `std::map`
сортирует ключи по первому элементу, а затем по второму. Это позволяет при итоговом
проходе по словарю получать элементы сразу в том порядке, который требует формат CCS
(сначала по столбцам, внутри столбца — по строкам).

### Вычислительное ядро

```cpp
void ProcessChunk(int start_idx, int end_idx, const CCSMatrix &mat_a, const std::vector<int> &a_cols,
                  const CCSMatrix &mat_b_t, LocalDict &local_map) {
  for (int k = start_idx; k < end_idx; ++k) {
    const int row_left = mat_a.row_index[static_cast<size_t>(k)];
    const int col_left = a_cols[static_cast<size_t>(k)];
    const std::complex<double> val_left = mat_a.values[static_cast<size_t>(k)];

    const int b_start = mat_b_t.col_index[static_cast<size_t>(col_left)];
    const int b_end = mat_b_t.col_index[static_cast<size_t>(col_left) + 1];

    for (int idx_p = b_start; idx_p < b_end; ++idx_p) {
      const int col_right = mat_b_t.row_index[static_cast<size_t>(idx_p)];
      const std::complex<double> val_right = mat_b_t.values[static_cast<size_t>(idx_p)];

      local_map[{col_right, row_left}] += val_left * val_right;
    }
  }
}
```

### Управление потоками

```cpp
const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;
  auto &mat_res = GetOutput();

  const int nnz_a = static_cast<int>(mat_a.values.size());
  if (nnz_a == 0) {
    mat_res.count_rows = mat_a.count_rows;
    mat_res.count_cols = mat_b.count_cols;
    mat_res.col_index.assign(static_cast<size_t>(mat_res.count_cols) + 1, 0);
    return true;
  }

  CCSMatrix mat_b_t;
  TransposeCCS(mat_b, mat_b_t);
  std::vector<int> a_cols = BuildColumnIndices(mat_a);

  unsigned int hw = std::thread::hardware_concurrency();
  int num_threads = std::min(nnz_a, (hw == 0) ? 1 : static_cast<int>(hw));

  std::vector<std::thread> threads;
  std::vector<LocalDict> local_maps(static_cast<size_t>(num_threads));
  int chunk = nnz_a / num_threads;
  int rem = nnz_a % num_threads;

  int start = 0;
  for (int i = 0; i < num_threads; ++i) {
    int end = start + chunk + (i < rem ? 1 : 0);
    threads.emplace_back(ProcessChunk, start, end, std::ref(mat_a), std::ref(a_cols), std::ref(mat_b_t),
                         std::ref(local_maps[static_cast<size_t>(i)]));
    start = end;
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  LocalDict global_map;
  AggregateLocalMaps(local_maps, global_map);
  ConstructResultMatrix(global_map, mat_a.count_rows, mat_b.count_cols, mat_res);

  return true;
```

## 5. Проверка корректности

Для верификации STL-версии использовались те же сценарии, что и в последовательной
реализации (единичные матрицы, случайные разреженные матрицы, прямоугольные матрицы).
Результаты полностью совпали с эталонными.

## 6. Экспериментальная среда

- **CPU:** Intel Core i5-12450H (8 ядер / 12 потоков)
- **RAM:** 16 GB
- **OS:** Ubuntu 24.04
- **Compiler:** GCC

## 7. Результаты производительности

Параметры теста: случайные разреженные матрицы 1000x1000 с плотностью заполнения 1%.

| Режим   | Потоки | Время Task_run (мс) | Ускорение | Эффективность |
| ------- | ------ | ------------------- | --------- | ------------- |
| SEQ     | 1      | 56.1                | 1.00x     | 100%          |
| STL     | 2      | 41.5                | 1.35x     | 67.5%         |
| STL     | 4      | 33.2                | 1.63x     | 40.7%         |

**Анализ результатов:**
Несмотря на параллельную работу, эффективность STL-версии ниже, чем у OMP и TBB.
Основная причина — интенсивное использование `std::map` и динамического выделения памяти
под узлы дерева внутри каждого потока. Объединение (слияние) нескольких словарей в конце
также является чисто последовательной и дорогой операцией.

## 8. Выводы

1. **Автономность:** Реализация полностью независима от внешних библиотек распараллеливания.
2. **Безопасность:** Стратегия "один поток — один словарь" позволила полностью избежать
   использования мьютексов и критических секций.
3. **Ограничения:** Использование ассоциативных контейнеров (`std::map`) внутри потоков
   создает значительный оверхед на аллокацию памяти, что ограничивает масштабируемость
   алгоритма на большом числе ядер.
4. **Упорядоченность:** Применение `std::pair<col, row>` в качестве ключа словаря
   значительно упростило фазу формирования финальной CCS-структуры.
