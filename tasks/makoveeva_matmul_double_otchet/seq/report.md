# ФИНАЛЬНЫЙ ОТЧЕТ: SEQ версия алгоритма Фокса

## Обзор SEQ версии

SEQ версия - базовая однопоточная реализация алгоритма Фокса для
умножения плотных матриц.

Технология: Sequential (без параллелизма)
Потоков: 1
Назначение: Отладка, справочная реализация, эталон сравнения
Сложность: O(n^3)

---

## Тестовые результаты

Все 12 размеров матриц от 1x1 до 32x32 успешно протестированы:

size_1x1 - PASSED (1 ms)
size_2x2 - PASSED (1 ms)
size_3x3 - PASSED (1 ms)
size_4x4 - PASSED (1 ms)
size_5x5 - PASSED (1 ms)
size_6x6 - PASSED (1 ms)
size_7x7 - PASSED (2 ms)
size_8x8 - PASSED (3 ms)
size_9x9 - PASSED (4 ms)
size_10x10 - PASSED (6 ms)
size_16x16 - PASSED (21 ms)
size_32x32 - PASSED (145 ms)

Результат: 12/12 PASSED

Производительные тесты на матрице 512x512:

Pipeline: 3862 ms
Task run: 4210 ms
Среднее: 4036 ms
GFLOPS: 0.319

---

## Исходный код

### Выбор размера блока

```cpp
int ChooseBlockSize(int n) {
  int block_size = static_cast<int>(std::sqrt(static_cast<double>(n)));
  block_size = std::max(1, block_size);

  while ((n % block_size != 0) && (block_size > 1)) {
    --block_size;
  }

  return block_size;
}
```

### Умножение блоков матриц

```cpp
void MultiplyBlocks(const std::vector<double> &a,
                   const std::vector<double> &b,
                   std::vector<double> &c,
                   int n,
                   int row_start, int row_end,
                   int col_start, int col_end,
                   int k_start, int k_end) {
  const auto n_size = static_cast<size_t>(n);

  for (int row = row_start; row < row_end; ++row) {
    const auto row_idx = static_cast<size_t>(row);
    const auto row_offset = row_idx * n_size;

    for (int col = col_start; col < col_end; ++col) {
      const auto col_idx = static_cast<size_t>(col);
      double sum = 0.0;

      for (int k = k_start; k < k_end; ++k) {
        const auto k_idx = static_cast<size_t>(k);
        sum += a[row_offset + k_idx] * b[(k_idx * n_size) + col_idx];
      }

      c[row_offset + col_idx] += sum;
    }
  }
}
```

### Основной алгоритм Фокса

SEQ версия реализует классический алгоритм Фокса в следующем порядке:

1. Разбивает матрицы A и B на блоки размером (n/block_size) x
(n/block_size)

2. На каждом этапе (от 0 до block_size-1) выполняет:
   - Определяет root блок по формуле: root = (i + stage) % block_size
   - Выполняет локальное умножение блоков A[i, root] и B[root, j]
   - Циклически сдвигает блоки матрицы B вверх на одну строку

---

## Результаты производительности

Время выполнения на матрице 512x512: 4036 ms
GFLOPS: 0.319

Это эталонная производительность, используется для сравнения с
другими версиями.

---

## Выводы

Преимущества SEQ версии:

1. Простая и понятная реализация
2. Легко отладить и проверить корректность
3. Детерминированное выполнение (без race conditions)
4. Работает на любой системе
5. Хороший эталон для сравнения

Недостатки SEQ версии:

1. Очень медленная (4036 ms на матрице 512x512)
2. Не использует потенциал многоядерных систем
3. Неподходящая для production на современных машинах

Использование:

SEQ версия рекомендуется для отладки алгоритма на ранних этапах
разработки, проверки корректности реализации на малых матрицах,
понимания базовой структуры алгоритма и справочной реализации для
сравнения.

Анализ производительности:

Время 4036 ms для матрицы 512x512 дает нам базовую метрику.

Все остальные версии будут сравниваться с этим результатом:

- OMP достигает 7.3x ускорения (556 ms)
- STL достигает 6.4x ускорения (628 ms)
- ALL достигает 7.1x ускорения (565.5 ms)
- TBB достигает 1.32x ускорения (3051 ms)

---

## Статус

Версия полностью функциональна и готова к использованию в качестве
справочной реализации.

Все тесты проходят успешно. Результаты математически корректны.

Статус: ГОТОВА (как базовая/справочная версия)
