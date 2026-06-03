# Построение выпуклой оболочки для компонент бинарного изображения - OpenMP

- **Student:** Дорогин Вадим Антонович, группа 3823Б1ПР3  
- **Technology:** OpenMP  
- **Вариант:** № 30

---

## 1. Контекст

Параллелизация переносится на этап **построения оболочек по уже найденным компонентам**:
поиск компонент (DFS) остаётся последовательным, как в SEQ; **`#pragma omp parallel for`**
обрабатывает индексы компонент независимо.

## 2. Постановка задачи

Постановка **совпадает** с корневым `report.md` и веткой **SEQ**; меняется
только реализация `RunImpl` (OpenMP на цикле по `components`).

**Вход / выход / проверка:** `BinaryImage`, **5** функциональных кейсов,
сравнение `convex_hulls` с эталоном.

## 3. Базовый алгоритм

1. `ThresholdImage` и `FindComponents` - последовательно (идентично SEQ).
2. `convex_hulls.resize(components.size())`.
3. Параллельный цикл: для каждой компоненты `BuildHull` или копия при `size <= 2`.

## 4. Схема распараллеливания

- **`num_threads(ppc::util::GetNumThreads())`** - число потоков из окружения курса.
- **`#pragma omp parallel for default(none) shared(components, convex_hulls)`** -
каждая итерация `i` пишет только в `convex_hulls[i]`, гонок нет.
- Неявный барьер в конце `parallel for` перед `GetOutput() = w_`.

Фрагмент:

```57:76:tasks/dorogin_v_bin_img_conv_hull/omp/src/ops_omp.cpp
bool DoroginVBinImgConvHullOMP::RunImpl() {
  FindComponents();
  w_.convex_hulls.clear();
  w_.convex_hulls.resize(w_.components.size());
#pragma omp parallel for default(none) shared(components, convex_hulls) num_threads(ppc::util::GetNumThreads())
  for (std::size_t i = 0; i < components.size(); ++i) {
    // BuildHull или копия компоненты
  }
  GetOutput() = w_;
  return true;
}
```

## 5. Детали реализации

- **Файлы:** `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`, класс `DoroginVBinImgConvHullOMP`.
- Подключены `<omp.h>` и `util/include/util.hpp` для `GetNumThreads()`.

## 6. Проверка корректности

Те же **5** кейсов, что и для SEQ; имя теста в gtest содержит `_omp_`.

## 7. Экспериментальная среда

| Параметр | Значение |
| -------- | -------- |
| Процессор | Intel Core i5-10400F, 6 ядер / 12 потоков, 2.90 ГГц |
| ОС | Microsoft Windows, 64-bit (сборка 10.0.26200) |
| ОЗУ | 16 ГБ DDR4, 2666 МТ/с |
| Устройство | DESKTOP-JQ9KQ17 |
| Сборка | **Release** |

- **`PPC_NUM_THREADS`** и **`OMP_NUM_THREADS`** - одинаковое значение (раннер курса).
- **`PPC_NUM_PROC=1`**.

```bash
export PPC_NUM_THREADS=4
export PPC_NUM_PROC=1
export OMP_NUM_THREADS=4
scripts/run_tests.py --running-type=threads
```

Perf: вход **600×600** (`kSize` в `tests/performance/main.cpp`).

## 8. Результаты

Заполнить после локальных прогонов `ppc_perf_tests` (режимы **`task_run`** / **`pipeline`**):

| Threads N | T_omp, с | S = T_seq / T_omp | Eff = S/N |
| --------- | -------- | ----------------- | --------- |
| 2         | 0,000860 | 0,71              | 0,36      |
| 4         | 0,000825 | 0,74              | 0,19      |
| 8         | 0,000948 | 0,65              | 0,08      |

`T_seq = 0,000614` с. Вход **600×600**, режим **`task_run`**, ПК DESKTOP-JQ9KQ17 (Release).

Сводные таблицы OMP/TBB/STL - корневой **`report.md`**, §5.

## 9. Выводы

OpenMP даёт простую параллелизацию независимых оболочек.
Узкое место при большом числе мелких компонент - накладные расходы OpenMP и
последовательный `FindComponents`. На изображении **600×600** с
диагоналями выигрыш ожидается при достаточном числе крупных компонент и `N > 1`.
Итоговое сравнение технологий - в корневом отчёте.

## 10. Источники

1. Документация курса.  
2. [OpenMP specification](https://www.openmp.org/).  
3. [cppreference - OpenMP](https://en.cppreference.com/w/cpp/keyword/omp_parallel).
