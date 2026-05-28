# Построение выпуклой оболочки для компонент бинарного изображения - oneTBB

- **Student:** Дорогин Вадим Антонович, группа 3823Б1ПР3  
- **Technology:** oneTBB  
- **Вариант:** № 30

---

## 1. Контекст

Используется **`tbb::parallel_for`** по индексам компонент для построения выпуклых оболочек.
Бинаризация и поиск компонент выполняются последовательно (как в SEQ/OMP).

## 2. Постановка задачи

Постановка **идентична** SEQ; в ветке TBB меняется параллельный этап в `RunImpl`.

## 3. Базовый алгоритм

Monotonic chain для оболочки; параллельный цикл по `components.size()`.

## 4. Схема распараллеливания

- `tbb::parallel_for` от `0` до `components.size()` с лямбдой по индексу `i`.
- Каждая итерация читает `components[i]` и пишет в `convex_hulls[i]`.
- Планировщик TBB распределяет итерации по потокам (`PPC_NUM_THREADS` в окружении).

Фрагмент:

```55:73:tasks/dorogin_v_bin_img_conv_hull/tbb/src/ops_tbb.cpp
  tbb::parallel_for(std::size_t{0}, components.size(), [&](std::size_t i) {
    const auto &comp = components[i];
    if (comp.empty()) {
      return;
    }
    if (comp.size() <= 2) {
      convex_hulls[i] = comp;
    } else {
      convex_hulls[i] = BuildHull(comp);
    }
  });
```

## 5. Детали реализации

- **Файлы:** `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`, класс `DoroginVBinImgConvHullTBB`.
- Заголовок: `oneapi/tbb/parallel_for.h`.

## 6. Проверка корректности

Общие **5** функциональных тестов; фильтр gtest `*tbb*`.

## 7. Экспериментальная среда

| Параметр | Значение |
| -------- | -------- |
| Процессор | Intel Core i5-10400F, 6 ядер / 12 потоков, 2.90 ГГц |
| ОС | Microsoft Windows, 64-bit (сборка 10.0.26200) |
| ОЗУ | 16 ГБ DDR4, 2666 МТ/с |
| Устройство | DESKTOP-JQ9KQ17 |
| Сборка | **Release** |

```bash
export PPC_NUM_THREADS=4
export PPC_NUM_PROC=1
export OMP_NUM_THREADS=4
scripts/run_tests.py --running-type=threads
```

## 8. Результаты

| N | T_tbb, с | S vs SEQ | Eff = S/N |
| - | -------- | -------- | --------- |
| 2 | 0,000688 | 0,89     | 0,45      |
| 4 | 0,000695 | 0,88     | 0,22      |
| 8 | 0,000861 | 0,71     | 0,09      |

`T_seq = 0,000614` с. Вход **600×600**, режим **`task_run`**.

Значения - из `ppc_perf_tests`, режим **`task_run`**, вход **600×600**. Сводка - корневой `report.md`.

## 9. Выводы

TBB удобен для параллельного цикла по компонентам без явного управления потоками.
По сравнению с OMP возможны отличия в планировании при неравномерных размерах компонент.
Абсолютные ускорения фиксируются локально и вносятся в корневой отчёт.

## 10. Источники

1. Документация курса.  
2. [oneAPI oneTBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html).  
3. [cppreference](https://en.cppreference.com/).
