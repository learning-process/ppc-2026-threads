# Построение выпуклой оболочки для компонент бинарного изображения - SEQ

- **Student:** Дорогин Вадим Антонович, группа 3823Б1ПР3  
- **Technology:** SEQ  
- **Вариант:** № 30

---

## 1. Контекст

Нужна **последовательная** реализация: бинаризация изображения, выделение
связных компонент (4-связность) и построение **выпуклой оболочки** для каждой компоненты.
Ветка SEQ - эталон корректности и базовое время для сравнения с OMP, TBB, STL и ALL.

## 2. Постановка задачи

**Цель.** По бинарному изображению найти все связные компоненты переднего плана и
для каждой построить выпуклую оболочку (множество вершин на границе выпуклости).

**Типы.** В `common/include/common.hpp`: **`InType` = `OutType` = `BinaryImage`**,
**`BaseTask` = `ppc::task::Task<InType, OutType>`**.
Структура содержит `width`, `height`, `pixels`, а также поля `components` и `convex_hulls`.

**Вход.** `BinaryImage` с `width > 0`, `height > 0` и `pixels.size() == width * height`.
Пиксели - значения яркости `uint8_t`.

**Выход.** Тот же `BinaryImage`, в `convex_hulls` - вектор оболочек
(каждая оболочка - `std::vector<Point>`).

**Проверка корректности.** Функциональные тесты в `tests/functional/main.cpp`
(**5** кейсов): сравнение `out.convex_hulls` с эталоном после нормализации вершин.

**Особые случаи.** Компонента из 1–2 точек возвращается без построения цепи;
пустые компоненты пропускаются.

## 3. Базовый алгоритм

1. **Порог:** `pixels[i] > 128` - 255, иначе 0 (`ThresholdImage`).
2. **Компоненты:** обход изображения; для каждого непосещённого белого пикселя - **DFS**
со стеком (`ExploreComponent`), 4 направления.
3. **Оболочка:** для компоненты с ≥ 3 точками - **monotonic chain**
после сортировки и удаления дубликатов (`BuildHull`); иначе копия компоненты.

**Сложность.** Пусть `N = width * height`, `K` - число компонент, `P_k` - размер k-й компоненты.
Бинаризация `O(N)`, поиск компонент `O(N)`, суммарно по оболочкам `O(Σ P_k log P_k)`.

## 4. Детали реализации

- **Файлы:** `seq/include/ops_seq.hpp`, `seq/src/ops_seq.cpp`,
класс `DoroginVBinImgConvHullSeq`.
- **`ValidationImpl`:** проверка размеров и длины `pixels`.
- **`PreProcessingImpl`:** копия входа и бинаризация.
- **`RunImpl`:** `FindComponents`, затем последовательное заполнение `convex_hulls`.
- **`PostProcessingImpl`:** без дополнительных действий; результат уже в `GetOutput()`.

Фрагмент построения оболочки:

```148:177:tasks/dorogin_v_bin_img_conv_hull/seq/src/ops_seq.cpp
std::vector<Point> DoroginVBinImgConvHullSeq::BuildHull(const std::vector<Point> &points) {
  std::vector<Point> pts = points;
  std::ranges::sort(pts, std::less<>{});
  const auto uniq = std::ranges::unique(pts);
  pts.erase(uniq.begin(), uniq.end());
  // ... нижняя и верхняя цепи, Cross ...
}
```

## 5. Проверка корректности

- **5** статических кейсов: одна точка, две компоненты, вертикальная линия, прямоугольник, «ромб» (манхэттен-окрестность).
- Сравнение через `Normalize` / `NormalizeAll` в тестах (сортировка вершин оболочки).

## 6. Экспериментальная среда

| Параметр | Значение |
| -------- | -------- |
| Процессор | Intel Core i5-10400F, 6 ядер / 12 потоков, 2.90 ГГц |
| ОС | Microsoft Windows, 64-bit (сборка 10.0.26200) |
| ОЗУ | 16 ГБ DDR4, 2666 МТ/с |
| Видеоадаптер | NVIDIA GeForce RTX 2060, 6 ГБ |
| Устройство | DESKTOP-JQ9KQ17 |
| Компилятор | согласно конфигурации CMake проекта PPC (Release) |
| Сборка | **Release** |

**Переменные:** `PPC_NUM_THREADS=1`, `PPC_NUM_PROC=1` для
согласованности с раннером.

**Сборка и func-тесты:**

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --parallel
export PPC_NUM_THREADS=1
export PPC_NUM_PROC=1
scripts/run_tests.py --running-type=threads
```

**Производительность:** `tests/performance/main.cpp`, изображение **600×600**,
диагонали и дополнительные точки в `SetUp`; режимы **`task_run`** и **`pipeline`**.

## 7. Результаты

### 7.1. Режим `task_run`

| Режим                | Время, с     | Роль             |
| -------------------- | ------------ | ---------------- |
| `task_run`, 600×600  | **0,000614** | Baseline `T_seq` |

Замер: `ppc_perf_tests`, фильтр `*task_run*dorogin_v_bin_img_conv_hull*seq*`,
`PPC_NUM_THREADS=1`.

### 7.2. Режим `pipeline`

В замер входит полный конвейер `Validation` - `PreProcessing` - `Run` - `PostProcessing`.
Сводка по всем backend-ам - в корневом `report.md`.

## 8. Выводы

SEQ задаёт корректную последовательную схему «бинаризация - компоненты - оболочки»
и служит эталоном для остальных технологий. На func-тестах (малые изображения)
время укладывается в лимит **`PPC_TASK_MAX_TIME` = 1 с** на один запуск пайплайна.
Для perf на **600×600** абсолютное время записывается локально и используется
в корневом отчёте для расчёта ускорений.

## 9. Источники

1. Документация курса «Параллельное программирование».  
2. Andrew, *Computational Geometry* (convex hull, monotonic chain).  
3. [cppreference](https://en.cppreference.com/) - `std::stack`, `std::ranges`.
