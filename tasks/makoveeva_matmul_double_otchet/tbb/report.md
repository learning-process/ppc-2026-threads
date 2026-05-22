
## TBB ВЕРСИЯ

### Обзор

TBB версия - реализация алгоритма Фокса с использованием Intel
Threading Building Blocks для параллельного выполнения.

Технология: Intel TBB
Потоков: 8
Ускорение: 1.32x (медленнее SEQ!)
Статус: Требует оптимизации

### Тестовые результаты

Все 12 тестов размеров от 1x1 до 32x32 прошли успешно.

Результат: 12/12 PASSED

Производительные тесты на матрице 512x512:

Pipeline: 2814 ms
Task run: 3289 ms
Среднее: 3051 ms
GFLOPS: 0.539

### Исходный код

```cpp
bool MatmulDoubleTBBTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const size_t n = n_;
  const auto &a = A_;
  const auto &b = B_;
  auto &c = C_;

  const size_t block_size = SelectBlockSize(n);

  if (n % block_size != 0) {
    return RunSimpleMultiply();
  }

  const size_t grid_size = n / block_size;

  oneapi::tbb::mutex write_mutex;

  oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, grid_size *
                                         grid_size * grid_size),
      [&](const oneapi::tbb::blocked_range<size_t> &range) {
        for (size_t step_i_j = range.begin();
             step_i_j != range.end(); ++step_i_j) {
          const size_t step = step_i_j / (grid_size * grid_size);
          const size_t i = (step_i_j % (grid_size * grid_size)) /
                          grid_size;
          const size_t j = step_i_j % grid_size;

          const size_t root = (i + step) % grid_size;

          std::vector<double> local_block(block_size * block_size,
                                         0.0);

          ComputeBlock(a, b, local_block, i, j, root,
                      block_size, n);

          AccumulateResult(c, local_block, i, j, block_size, n,
                          write_mutex);
        }
      });

  GetOutput() = C_;
  return true;
}
```

### Выводы

TBB версия функционирует корректно, но производительность неоптимальна.

Причины медленности:

- Overhead синхронизации (mutex используется часто)
- Выделение памяти в каждой итерации
- Слишком мелкая гранулярность задач
- Task scheduler overhead

Рекомендуется использовать для более сложных паттернов параллелизма
после оптимизации реализации.

Статус: ГОТОВА, требует доработки для production

---
