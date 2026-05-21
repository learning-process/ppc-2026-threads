# Построение выпуклой оболочки для компонент бинарного изображения - ALL (MPI + TBB)

- **Student:** Дорогин Вадим Антонович, группа 3823Б1ПР3  
- **Technology:** ALL (MPI + TBB)  
- **Вариант:** № 30

---

## 1. Контекст

**Гибридная** реализация: вычисления на **rank 0** (бинаризация, компоненты, оболочки),
затем рассылка результата всем рангам через **MPI**.
Внутри rank 0: **TBB** на бинаризации и на построении оболочек; поиск компонент - **последовательный BFS** (очередь).

## 2. Постановка задачи

**Цель.** Та же, что в SEQ: для бинарного изображения получить `convex_hulls`
для всех связных компонент.

**Особенность ALL.** После вычисления на rank 0 результат **сериализуется** (`FlattenHulls`)
и передаётся **`MPI_Bcast`**; все ранги восстанавливают `convex_hulls` (`RestoreHulls`).
 Флаг валидности также рассылается (`MPI_Bcast` для `ok`), чтобы избежать зависания при ошибке
 на rank 0.

**Запуск тестов.** Func/perf для `_all_` выполняются под **`mpirun`**
(`scripts/run_tests.py --running-type=processes`); без MPI кейсы пропускаются.

## 3. Базовый алгоритм

1. `ThresholdPixelsParallel` - TBB по `pixels`.
2. `CollectComponentsSequential` - BFS, 4-связность.
3. `tbb::parallel_for` по компонентам - `BuildHull` или копия.
4. `FlattenHulls` - `MPI_Bcast` размера буфера и данных - `RestoreHulls` на всех рангах.
5. `MPI_Barrier`, `GetOutput() = work_`.

## 4. Межпроцессная схема

| Шаг | Операция |
| --- | -------- |
| 1 | `MPI_Comm_rank` |
| 2 | На rank 0: `ValidationImpl` - `ok`, `MPI_Bcast(&ok, …)` |
| 3 | Вычисление на rank 0 |
| 4 | `MPI_Bcast` длины и массива `packed` |
| 5 | `RestoreHulls`, очистка `components` |

## 5. Внутрипроцессная схема

- TBB: бинаризация и оболочки (как в TBB-ветке по смыслу).
- BFS вместо DFS (эквивалентная 4-связность для компонент).

## 6. Детали реализации

- **Файлы:** `all/include/ops_all.hpp`, `all/src/ops_all.cpp`, класс `DoroginVBinImgConvHullALL`.
- Сериализация: число оболочек, для каждой - число точек и пары `(x, y)`.

Фрагмент MPI:

```197:240:tasks/dorogin_v_bin_img_conv_hull/all/src/ops_all.cpp
bool DoroginVBinImgConvHullALL::RunImpl() {
  // rank 0: CollectComponentsSequential + tbb::parallel_for по hulls
  // FlattenHulls, MPI_Bcast packed_len и packed, RestoreHulls
  MPI_Barrier(MPI_COMM_WORLD);
  GetOutput() = work_;
  return true;
}
```

## 7. Проверка корректности

- **5** функциональных кейсов × ALL (под MPI).
- Сравнение `convex_hulls` с эталоном на каждом ранге после `RestoreHulls`.

## 8. Экспериментальная среда

| Параметр | Значение |
| -------- | -------- |
| Процессор | Intel Core i5-10400F, 6 ядер / 12 потоков, 2.90 ГГц |
| ОС | Microsoft Windows, 64-bit (сборка 10.0.26200) |
| ОЗУ | 16 ГБ DDR4, 2666 МТ/с |
| Устройство | DESKTOP-JQ9KQ17 |
| Сборка | **Release** |

```bash
export PPC_NUM_THREADS=4
export PPC_NUM_PROC=2
export OMP_NUM_THREADS=4
scripts/run_tests.py --running-type=processes
```

Пример точечного perf (уточнить имя через `--gtest_list_tests`):

```bash
mpirun -np 1 ./build/bin/ppc_perf_tests \
  '--gtest_filter=*dorogin*all*task_run*' --gtest_brief=1
```

## 9. Результаты

### 9.1. `P = 1` (один MPI-процесс)

| N потоков | T_all, с | S vs SEQ | Eff = S/N |
| --------- | -------- | -------- | --------- |
| 2         | 0,000815 | 0,75     | 0,38      |
| 4         | 0,000807 | 0,76     | 0,19      |
| 8         | 0,000953 | 0,64     | 0,08      |

`T_seq = 0,000614` с, `mpiexec -n 1`.

### 9.2. `P > 1`

| P | N | T_all, с | S vs SEQ | Eff = S/(P·N) |
| - | - | -------- | -------- | ------------- |
| 2 | 2 | 0,001071 | 0,57     | 0,14          |
| 2 | 4 | 0,000906 | 0,68     | 0,08          |

`mpiexec -n 2`, `T_seq = 0,000614` с.

Сводка - корневой **`report.md`**.

## 10. Выводы

ALL объединяет TBB на rank 0 и MPI для репликации результата. При **`P = 1`** измеряется
в основном TBB + накладные расходы MPI; при **`P > 1`** на одном узле доминируют коммуникация
и конкуренция за ядра. Для учебной постановки на одной машине разумно
сравнивать с OMP/TBB/STL при **`P = 1`**.

## 11. Источники

1. Документация курса.  
2. [MPI Forum](https://www.mpi-forum.org/).  
3. [oneTBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html).
