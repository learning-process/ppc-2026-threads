# Обработка контуров бинарных компонент — ALL (MPI + OpenMP)

- **Student:** Перяшкин Василий Андреевич
- **Technology:** ALL (MPI + OpenMP)
- **Variant:** 30
- **Group:** 3823Б1ПР4

## 1. Контекст

Гибридная версия комбинирует два уровня параллелизма:

- **Межпроцессный (MPI):** компоненты распределяются между рангами через `MPI_Scatterv`;
  результаты собираются через `MPI_Gatherv` и рассылаются через `MPI_Bcast`.
- **Внутрипроцессный (OpenMP):** каждый ранг строит оболочки параллельно
  через `#pragma omp parallel for`.

BFS выполняется только на ранге 0. Подход оправдан, когда число компонент велико и
стоимость передачи данных окупается параллельным вычислением оболочек.

## 2. Постановка задачи

Вход, выход и ограничения идентичны SEQ. Все ранги завершают `SolveALL` с одинаковым
набором оболочек (благодаря финальному `MPI_Bcast`).
Число MPI-процессов: `PPC_NUM_PROC`; потоков на ранг: `PPC_NUM_THREADS`.

## 3. Базовый алгоритм

BFS + монотонная цепь Эндрю.

## 4. Межпроцессная схема

Весь конвейер `SolveALL` разделён на 5 этапов:

```text
Ранг 0                          Ранги 1..P-1
──────────────────              ──────────────────
1. ExtractComponents4           (ждут scatter)
2. FlattenDistributedComponents
   MPI_Scatter(send_counts)  →→ recv_count
   MPI_Scatterv(flat_send)   →→ flat_recv
                                3. UnflattenComponents
                                   BuildLocalHulls (OMP)
                                   FlattenComponents
   MPI_Gatherv             ←←  flat_local_hulls
4. UnflattenComponents (root)
   MPI_Bcast(broadcast_size) →→ broadcast_size
   MPI_Bcast(flat_broadcast) →→ flat_broadcast
                                5. UnflattenComponents
```

**Роли рангов:**

| Ранг 0                             | Ранги 1..P-1                         |
|------------------------------------|--------------------------------------|
| BFS (извлечение компонент)         | Ждут `MPI_Scatter`                   |
| Сериализация в плоский массив      | Получают плоский массив компонент    |
| Десериализуют и строят оболочки    | Строят оболочки локально (OMP)       |
| Собирает результаты (`Gatherv`)    | Отправляют результаты                |
| Рассылает итог (`Bcast`)           | Получают итог                        |

**Сериализация компонент:** `FlattenComponents` упаковывает компоненты в плоский `int`-массив:
`[size_i, x_0, y_0, ..., x_{n-1}, y_{n-1}]`. `UnflattenComponents` обращает процедуру.

**Балансировка:** ранг `i` получает `total/P + (i < total%P ? 1 : 0)` компонент.

**MPI-коммуникации:**

| Вызов                        | Назначение                                                    |
|------------------------------|---------------------------------------------------------------|
| `MPI_Scatter(send_counts)`   | Сообщить каждому рангу размер его порции плоского массива     |
| `MPI_Scatterv(flat_send)`    | Разослать плоские данные компонент с переменными смещениями   |
| `MPI_Gather(local_size)`     | Собрать размеры результатов у ранга 0                         |
| `MPI_Gatherv(flat_hulls)`    | Собрать плоские оболочки с переменными смещениями             |
| `MPI_Bcast(size + data)`     | Разослать итоговый результат всем рангам                      |

## 5. Внутрипроцессная схема

Внутри каждого ранга применяется `BuildLocalHulls`:

```cpp
// all/src/ops_all.cpp  (строки 210–221)
inline OutType BuildLocalHulls(std::vector<std::vector<Point>> local_components) {
  OutType local_hulls(local_components.size());
  const int n = static_cast<int>(local_components.size());

#pragma omp parallel for default(none) shared(local_components, local_hulls, n)
  for (int i = 0; i < n; ++i) {
    const auto idx = static_cast<std::size_t>(i);
    local_hulls[idx] = ConvexHullMonotonicChain(std::move(local_components[idx]));
  }
  return local_hulls;
}
```

`omp parallel for` без явного `schedule` — по умолчанию `schedule(static)`.
Каждый ранг сериализует `local_hulls` через `FlattenComponents` перед `MPI_Gatherv`.

## 6. Детали реализации

**Файлы:** [all/include/ops_all.hpp](../all/include/ops_all.hpp),
[all/src/ops_all.cpp](../all/src/ops_all.cpp)

**Потенциальные узкие места:**

1. `MPI_Bcast` финального результата — при большом числе компонент занимает ощутимое время.
2. Сериализация/десериализация (`FlattenComponents`/`UnflattenComponents`) — последовательна,
   вносит overhead при большом числе компонент.
3. BFS на ранге 0 — последовательный узкий узел как и во всех других версиях.

**Точки синхронизации между уровнями:**

- После `MPI_Scatterv` — неявная синхронизация: ранги не начинают вычисление оболочек,
  пока не получат данные.
- После `MPI_Gatherv` — неявная синхронизация: ранг 0 не рассылает результат до получения.
- Явных `MPI_Barrier` нет — коллективные операции сами являются точками синхронизации.

## 7. Проверка корректности

Все ранги завершают с идентичным `local_out_` благодаря `MPI_Bcast`.
Функциональные тесты проверяют совпадение с SEQ при запуске с 2 и 4 рангами:

```bash
PPC_NUM_PROC=2 PPC_NUM_THREADS=2 scripts/run_tests.py --running-type=processes \
  --gtest_filter="*peryashkin_v*all*"
```

Согласованность между рангами гарантируется схемой: все ранги получают одни данные
через `MPI_Bcast`, оболочки вычисляются детерминированно.

## 8. Экспериментальная среда

| Параметр         | Значение                                      |
|------------------|-----------------------------------------------|
| CPU              | Apple M2 (4P + 4E cores, ARM64)               |
| Ядра / потоки    | 8 / 8                                         |
| RAM              | 8 GB                                          |
| OS               | macOS 13″ / Docker (Linux container)          |
| Компилятор       | GCC 13.3.0, `-std=c++20`, `-O3`               |
| MPI              | Open MPI 4.1.6                                |
| Build type       | Release (`-O3 -DNDEBUG`)                      |
| PPC_NUM_PROC     | 2, 4                                          |
| PPC_NUM_THREADS  | 1, 2 (на ранг)                                |

**Команды запуска:**

```bash
# 2 ранга × 2 потока на ранг = 4 worker'а
PPC_NUM_PROC=2 PPC_NUM_THREADS=2 scripts/run_tests.py --running-type=processes \
  --gtest_filter="*peryashkin_v*all*"

# Производительность
PPC_NUM_PROC=4 PPC_NUM_THREADS=1 scripts/run_tests.py --running-type=performance \
  --gtest_filter="*peryashkin_v*all*"
```

## 9. Результаты

| Ranks | Threads/rank | Workers | Медианное время (task) | Медианное время (pipeline) | Speedup task | Speedup pipe |
|-------|--------------|---------|------------------------|----------------------------|--------------|--------------|
| 2     | 1            | 2       | 1.379 мс               | 1.521 мс                   | 0.99×        | 0.95×        |
| 4     | 1            | 4       | 2.565 мс               | 2.639 мс                   | 0.53×        | 0.55×        |
| 2     | 2            | 4       | 1.393 мс               | 1.472 мс                   | 0.98×        | 0.99×        |

> Эффективность считается по числу workers = ranks × threads_per_rank.
>
> **Ожидаемое поведение:** при малом числе компонент гибридная версия может оказаться
> медленнее SEQ из-за сериализации и MPI-коммуникаций. Выигрыш — при больших изображениях
> с сотнями тысяч компонент.

## 10. Выводы

Гибридная схема оправдана при большом числе компонент и нескольких MPI-узлах.
На одной машине чистый OMP/TBB эффективнее: нет overhead на сериализацию и MPI-коммуникации.
Узкое место — BFS на ранге 0 и финальный `MPI_Bcast`.
