# Отчет STL: сортировка Хоара с простым слиянием

## Контекст и базовый алгоритм

STL-версия использует `std::thread`: массив делится на блоки по 64 элемента, каждый блок сортируется quicksort Хоара,
после чего уровни простого слияния выполняются параллельно по независимым парам диапазонов
([`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L16), [`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L143)). Вход и выход имеют тип
`std::vector<int>` ([`common/include/common.hpp`](../common/include/common.hpp#L11)).

## Разбиение, буферы и синхронизация

Количество потоков выбирается внутри реализации как `min(task_count, hardware_concurrency)`; переменная окружения
`PPC_NUM_THREADS` в STL-коде не читается. Если `hardware_concurrency()` вернул `0`, используется fallback `2`
([`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L18)). На тестовой машине `nproc` вернул `12`, поэтому в таблице конфигурация
обозначена как `auto (12)`, а не как вручную заданные 12 потоков. `RunInThreads` распределяет задачи циклически:
`task_index += thread_count` ([`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L41)). Локальные буферы представлены стековыми
переменными lambda, общий `merged_data` безопасен, потому что каждый task пишет только свой диапазон
([`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L163)). `atomic` и `mutex` не используются, потому что общей изменяемой
скалярной структуры нет.

`join` вызывается после запуска всех `std::thread`, чтобы сначала создать весь пул рабочих потоков, а затем дождаться
завершения каждого. Если делать `join` сразу после `emplace_back`, выполнение стало бы последовательным по потокам.
Семантика `join` подтверждена cppreference: завершение потока синхронизируется с успешным возвратом из `join`
([cppreference `std::thread::join`](https://en.cppreference.com/w/cpp/thread/thread/join.html)).

Фрагмент, [`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L38): создание всех потоков и последующий join.

```cpp
std::vector<std::thread> threads;
threads.reserve(thread_count);

for (size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
  threads.emplace_back([thread_index, thread_count, task_count, &function]() {
    for (size_t task_index = thread_index; task_index < task_count; task_index += thread_count) {
      function(task_index);
    }
  });
}

for (auto &thread : threads) {
  thread.join();
}
```

## Детали pipeline

`ValidationImpl` проверяет непустой вход ([`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L133)). `PreProcessingImpl` копирует
вход в `data_` ([`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L137)). `RunImpl` сортирует блоки и сливает уровни
([`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L143)). `PostProcessingImpl` проверяет сортировку и записывает результат в
`GetOutput()` только при успехе ([`stl/src/ops_stl.cpp`](src/ops_stl.cpp#L182)).

## Корректность и среда

В исходном тесте STL-backend добавлен в общий список задач ([`tests/functional/main.cpp`](../tests/functional/main.cpp#L84)),
а эталон строится через `std::ranges::sort` ([`tests/functional/main.cpp`](../tests/functional/main.cpp#L35)). Свежий
запуск `build_olesnitskiy/bin/ppc_func_tests` подтвердил `seq/omp/stl/tbb`: 60 passed; ALL отдельно прошел под
`mpirun -np 2`: 15 passed.

## Результаты

Baseline: `seq` `TaskRun = 0.0058254364 s`; для pipeline baseline `0.0068995056 s`. Performance-вход описан в
[`tests/performance/main.cpp`](../tests/performance/main.cpp#L20). Framework выполняет 5 повторов по умолчанию
([`modules/performance/include/performance.hpp`](../../../modules/performance/include/performance.hpp#L21)).

| workers | time | speedup | efficiency | notes |
|---|---:|---:|---:|---|
| auto (12 на тестовой машине) | 0.0047718214 s | 1.221 | 0.102 | `TaskRun`; число workers взято из `hardware_concurrency()` |
| auto (12 на тестовой машине) | 0.0055675704 s | 1.239 | 0.103 | `pipeline`; запуск с `PPC_NUM_THREADS=1`, но STL-код его не использует |
| auto (12 на тестовой машине) | 0.0081357990 s | 0.848 | 0.071 | `pipeline`; запуск с `PPC_NUM_THREADS=4`, но STL-код его не использует |

## Выводы

По коду STL-версия имеет те же фазы, что TBB: независимая сортировка блоков и независимое слияние. На `N=100000`
получено `0.0047718214 s`, speedup `1.221`; efficiency низкая (`0.102`), потому что реализация автоматически создает
до `hardware_concurrency()` рабочих потоков, а объем работы мал для такого числа workers.
