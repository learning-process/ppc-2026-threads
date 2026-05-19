# Обработка контуров бинарных компонент — ALL (MPI + OpenMP)

- **Student:** Перяшкин Василий Андреевич
- **Technology:** ALL (MPI + OpenMP)
- **Variant:** 4
- **Group:** 3823Б1ПР4

## 1. Контекст

Гибридная версия комбинирует два уровня параллелизма:
- **Межпроцессный (MPI):** компоненты распределяются между MPI-рангами через `MPI_Scatterv`, результаты собираются через `MPI_Gatherv` и рассылаются через `MPI_Bcast`.
- **Внутрипроцессный (OpenMP):** каждый ранг вычисляет выпуклые оболочки своих компонент параллельно через `#pragma omp parallel for`.

BFS по-прежнему выполняется только на ранге 0. Данный подход оправдан, когда общее число компонент велико и стоимость сериализации/передачи данных окупается параллельным вычислением оболочек на нескольких узлах.

Базовый SEQ-эталон — [seq/report.md](../seq/report.md).

## 2. Постановка задачи

Вход, выход и ограничения идентичны SEQ. Все ранги завершают `SolveALL` с одним и тем же набором оболочек (благодаря финальному `MPI_Bcast`). Число MPI-процессов задаётся `PPC_NUM_PROC`, число OpenMP-потоков на ранг — `PPC_NUM_THREADS`.

## 3. Базовый алгоритм

BFS + монотонная цепь Эндрю. Подробнее — в [seq/report.md](../seq/report.md).

## 4. Межпроцессная схема

Весь конвейер `SolveALL` разделён на 5 этапов:

```
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

**Сериализация компонент:** структура `FlattenComponents` упаковывает вектор компонент в плоский `int`-массив: `[size_i, x_0, y_0, ..., x_{n-1}, y_{n-1}]` для каждой компоненты. `UnflattenComponents` обращает процедуру.

**Балансировка:** `MakeComponentCounts` распределяет компоненты равномерно: ранг `i` получает `total/P + (i < total%P ? 1 : 0)` компонент. Round-robin по равному числу компонент.

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
  const int local_components_count = static_cast<int>(local_components.size());

#pragma omp parallel for default(none) shared(local_components, local_hulls, local_components_count)
  for (int i = 0; i < local_components_count; ++i) {
    const auto idx = static_cast<std::size_t>(i);
    local_hulls[idx] = ConvexHullMonotonicChain(std::move(local_components[idx]));
  }
  return local_hulls;
}
```

Это то же самое `omp parallel for` без `schedule(static)`, то есть используется `schedule(static)` по умолчанию. Объединение локальных результатов: каждый ранг сериализует `local_hulls` через `FlattenComponents` перед `MPI_Gatherv`.

## 6. Детали реализации

**Файлы:** [all/include/ops_all.hpp](../all/include/ops_all.hpp), [all/src/ops_all.cpp](../all/src/ops_all.cpp)

**Потенциальные узкие места:**
1. `MPI_Bcast` финального результата рассылает все оболочки всем рангам — при большом числе компонент это может занимать ощутимое время.
2. Сериализация/десериализация (`FlattenComponents` / `UnflattenComponents`) выполняется последовательно — при очень большом числе компонент это тоже вносит overhead.
3. BFS на ранге 0 — последовательный узкий узел как и во всех других версиях.

**Точки синхронизации между уровнями:**
- После `MPI_Scatterv` — неявная синхронизация: ранги не начинают вычисление оболочек, пока не получат данные.
- После `MPI_Gatherv` — неявная синхронизация: ранг 0 не рассылает результат до получения от всех.
- Явных `MPI_Barrier` в коде нет — все коллективные операции сами являются точками синхронизации.

## 7. Проверка корректности

Все ранги завершают с идентичным `local_out_` благодаря финальному `MPI_Bcast`. Функциональные тесты проверяют совпадение с SEQ-эталоном при запуске с 2 и 4 рангами:
```bash
PPC_NUM_PROC=2 PPC_NUM_THREADS=2 scripts/run_tests.py --running-type=processes \
  --gtest_filter="*peryashkin_v*all*"
```

Согласованность результатов между рангами гарантируется самой схемой: все ранги получают одни и те же данные через `MPI_Bcast`, а оболочки вычисляются детерминированно.

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

SEQ baseline: task **1.370 мс**, pipeline **1.452 мс** (медиана 5 повторов, 512×512 checkerboard).

| Ranks | Threads/rank | Total workers | Медианное время (task) | Медианное время (pipeline) | Speedup task | Speedup pipe | Efficiency task | Efficiency pipe |
|-------|--------------|---------------|------------------------|----------------------------|--------------|--------------|-----------------|-----------------|
| 2     | 1            | 2             | 1.379 мс               | 1.521 мс                   | 0.99×        | 0.95×        | 50%             | 48%             |
| 4     | 1            | 4             | 2.565 мс               | 2.639 мс                   | 0.53×        | 0.55×        | 13%             | 14%             |
| 2     | 2            | 4             | 1.393 мс               | 1.472 мс                   | 0.98×        | 0.99×        | 25%             | 25%             |

> Эффективность считается по общему числу workers = ranks × threads_per_rank.

> **Ожидаемое поведение:** при малом числе компонент (или маленьком изображении) гибридная версия может оказаться медленнее SEQ из-за стоимости сериализации и MPI-коммуникаций. Выигрыш проявляется при больших изображениях с сотнями/тысячами компонент, когда стоимость передачи данных ниже выигрыша от параллельного вычисления оболочек.

## 10. Выводы

Гибридная схема оправдана при большом числе компонент и нескольких MPI-узлах. На одной машине с разделяемой памятью чистый OMP/TBB эффективнее из-за отсутствия overhead на сериализацию и MPI-коммуникации. Основное узкое место — последовательный BFS на ранге 0 и финальный `MPI_Bcast` результата.
