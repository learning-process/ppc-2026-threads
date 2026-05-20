# Сортировка Хоара с простым слиянием — STL

- **Студент:** Олесницкий Владимир Тарасович, 3823Б1ПР2
- **Технология:** STL
- **Вариант:** 13

## 1. Контекст

STL-версия использует `std::thread` для параллельной сортировки блоков и
последующего слияния независимых пар диапазонов. В отличие от OMP и TBB, число
рабочих потоков выбирается внутри реализации через
`std::thread::hardware_concurrency`.

## 2. Постановка задачи

- **Входные данные:** непустой объект `std::vector<int>`.
- **Выходные данные:** отсортированный по неубыванию `std::vector<int>`.
- **Baseline:** последовательная версия с временем `T_seq = 0.0058254364 s`.

## 3. Базовый алгоритм

Массив делится на блоки по 64 элемента. Каждый блок сортируется quicksort Хоара,
после чего уровни простого слияния выполняются параллельно по независимым парам
диапазонов.

## 4. Схема распараллеливания

Количество потоков выбирается как `min(task_count, hardware_concurrency)`. Если
`hardware_concurrency()` вернул `0`, используется fallback `2`. На тестовой
машине `nproc` вернул `12`, поэтому в результатах конфигурация обозначена как
`auto (12)`.

`RunInThreads` распределяет задачи циклически: `task_index += thread_count`.
Локальные буферы представлены стековыми переменными lambda, общий
`merged_data` безопасен, потому что каждая task пишет только свой диапазон.
`atomic` и `mutex` не используются, так как общей изменяемой скалярной
структуры нет.

Фрагмент создания и ожидания потоков:

```cpp
std::vector<std::thread> threads;
threads.reserve(thread_count);

for (size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
  threads.emplace_back([thread_index, thread_count, task_count, &function]() {
    for (size_t task_index = thread_index; task_index < task_count;
         task_index += thread_count) {
      function(task_index);
    }
  });
}

for (auto &thread : threads) {
  thread.join();
}
```

`join` вызывается после запуска всех `std::thread`, чтобы сначала создать весь
пул рабочих потоков, а затем дождаться завершения каждого. Если делать `join`
сразу после `emplace_back`, выполнение стало бы последовательным по потокам.

## 5. Детали реализации

`ValidationImpl` проверяет непустой вход. `PreProcessingImpl` копирует вход в
`data_`. `RunImpl` сортирует блоки и сливает уровни. `PostProcessingImpl`
проверяет сортировку и записывает результат в `GetOutput()` только при успехе.

## 6. Проверка корректности

STL-backend добавлен в общий список задач функционального теста, а эталон
строится через `std::ranges::sort`. Свежий запуск
`build_olesnitskiy/bin/ppc_func_tests` подтвердил `seq/omp/stl/tbb`: 60 passed.
ALL отдельно прошел под `mpirun -np 2`: 15 passed.

## 7. Экспериментальная среда

- **Сборка:** `build_olesnitskiy`
- **Compiler:** `g++-14`
- **Flags:** `-O3 -DNDEBUG`, `std=gnu++23`
- **Размер входных данных:** `N=100000`
- **Baseline TaskRun:** `0.0058254364 s`
- **Baseline pipeline:** `0.0068995056 s`
- **Число повторов:** 5 по умолчанию

Переменная `PPC_NUM_THREADS` в STL-коде не читается.

## 8. Результаты

- workers: auto (12 на тестовой машине); time: 0.0047718214 s; speedup: 1.221;
  efficiency: 0.102; notes: `TaskRun`; `hardware_concurrency()`.
- workers: auto (12 на тестовой машине); time: 0.0055675704 s; speedup: 1.239;
  efficiency: 0.103; notes: `pipeline`; env=1 не влияет.
- workers: auto (12 на тестовой машине); time: 0.0081357990 s; speedup: 0.848;
  efficiency: 0.071; notes: `pipeline`; env=4 не влияет.

## 9. Выводы

STL-версия имеет те же фазы, что TBB: независимая сортировка блоков и
независимое слияние. На `N=100000` получено `0.0047718214 s`, speedup `1.221`.
Efficiency низкая (`0.102`), потому что реализация автоматически создает до
`hardware_concurrency()` рабочих потоков, а объем работы мал для такого числа
workers.
