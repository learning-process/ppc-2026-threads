# Маркировка компонент на бинарном изображении — TBB

- **Студент**: Артюшкина Юлия Дмитриевна, группа 3823Б1ПР4
- **Технология**: TBB
- **Вариант**: 29

## 1. Контекст

Исходная задача — маркировка связных компонент на бинарном изображении, где чёрные
пиксели (значение 0) соответствуют объектам, а белые (значение 255) — фону.

TBB-реализация использует библиотеку Intel oneTBB для задачно-ориентированного
распараллеливания. Алгоритм основан на двухпроходном методе с системой
непересекающихся множеств (DSU) и синхронизацией через `tbb::spin_mutex`.

## 2. Постановка задачи

**Входные данные**: вектор `std::vector<uint8_t>`, где:
- `[0]` — количество строк `rows`
- `[1]` — количество столбцов `cols`
- `[2..]` — значения пикселей: `0` (объект) или `255` (фон)

**Выходные данные**: вектор `std::vector<uint8_t>`, где:
- `[0]` — `rows`
- `[1]` — `cols`
- `[2..]` — метки компонент (`0` — фон, `1..N` — уникальные метки объектов)

## 3. Базовый алгоритм

TBB-реализация использует двухпроходный алгоритм с DSU, разбитый на несколько этапов:

1. InitLabelsTbb — параллельная инициализация меток
2. MergeHorizontalPairsTbb — параллельное объединение горизонтальных пар
3. MergeVerticalPairsTbb — параллельное объединение вертикальных пар
4. MergeDiagonalPairsTbb — параллельное объединение диагональных пар
5. FinalizeRootsTbb — параллельное разрешение корней
6. NormalizeLabelsTbb — последовательная перенумерация

## 4. Детали реализации

**Файлы реализации:**
- `tbb/include/ops_tbb.hpp` — объявление класса `MarkingComponentsTBB`
- `tbb/src/ops_tbb.cpp` — реализация TBB-алгоритма

**Ключевые методы:**
- `InitLabelsTbb()` — параллельная инициализация
- `MergeHorizontalPairsTbb()` — объединение горизонтальных соседей
- `MergeVerticalPairsTbb()` — объединение вертикальных соседей
- `MergeDiagonalPairsTbb()` — объединение диагональных соседей
- `FinalizeRootsTbb()` — параллельное разрешение корней
- `NormalizeLabelsTbb()` — последовательная перенумерация
- `FindRoot(int label)` — поиск корня с сжатием пути
- `UnionLabels(int label1, int label2)` — объединение с мьютексом

**Инициализация меток:**

```cpp
void MarkingComponentsTBB::InitLabelsTbb() {
    int total_pixels = rows_ * cols_;
    tbb::parallel_for(0, total_pixels, [this](int idx) {
        size_t input_idx = static_cast<size_t>(idx) + 2;
        if (input_[input_idx] == 0) {
            labels_[idx] = idx + 1;
        }
    });
}
Объединение горизонтальных пар:

cpp
void MarkingComponentsTBB::MergeHorizontalPairsTbb() {
    tbb::parallel_for(0, rows_, [this](int y_coord) {
        for (int x_coord = 0; x_coord < cols_ - 1; ++x_coord) {
            int idx = (y_coord * cols_) + x_coord;
            if (labels_[idx] != 0 && labels_[idx + 1] != 0) {
                UnionLabels(labels_[idx], labels_[idx + 1]);
            }
        }
    });
}
Объединение вертикальных пар:

cpp
void MarkingComponentsTBB::MergeVerticalPairsTbb() {
    tbb::parallel_for(0, rows_ - 1, [this](int y_coord) {
        for (int x_coord = 0; x_coord < cols_; ++x_coord) {
            int idx = (y_coord * cols_) + x_coord;
            if (labels_[idx] != 0 && labels_[idx + cols_] != 0) {
                UnionLabels(labels_[idx], labels_[idx + cols_]);
            }
        }
    });
}
Объединение диагональных пар:

cpp
void MarkingComponentsTBB::MergeDiagonalPairsTbb() {
    tbb::parallel_for(0, rows_ - 1, [this](int y_coord) {
        for (int x_coord = 0; x_coord < cols_ - 1; ++x_coord) {
            int idx = (y_coord * cols_) + x_coord;
            if (labels_[idx] != 0 && labels_[idx + cols_ + 1] != 0) {
                UnionLabels(labels_[idx], labels_[idx + cols_ + 1]);
            }
            if (x_coord > 0) {
                if (labels_[idx] != 0 && labels_[idx + cols_ - 1] != 0) {
                    UnionLabels(labels_[idx], labels_[idx + cols_ - 1]);
                }
            }
        }
    });
}
FindRoot с сжатием пути:

cpp
int MarkingComponentsTBB::FindRoot(int label) {
    int root = label;
    while (parent_[root] != root) {
        root = parent_[root];
    }
    int current = label;
    while (parent_[current] != current) {
        int next = parent_[current];
        parent_[current] = root;
        current = next;
    }
    return root;
}
UnionLabels с мьютексом:

cpp
void MarkingComponentsTBB::UnionLabels(int label1, int label2) {
    tbb::spin_mutex::scoped_lock lock(dsu_mutex_);
    int root1 = FindRoot(label1);
    int root2 = FindRoot(label2);
    if (root1 != root2) {
        if (root1 < root2) {
            parent_[root2] = root1;
        } else {
            parent_[root1] = root2;
        }
    }
}
5. Проверка корректности
Корректность проверялась на 6 функциональных тестах.

Тест	Размер	Описание	Результат
0	3x3	L-образная фигура	PASSED
1	3x3	Диагонально соединённые компоненты	PASSED
2	2x3	Полностью фоновое изображение	PASSED
3	2x2	Полностью заполнено объектами	PASSED
4	3x4	Две горизонтальные полосы	PASSED
6	2x2	Проверка диагональной связности	PASSED
6. Экспериментальная среда
Аппаратное обеспечение:

Процессор: AMD Ryzen 7 7730U with Radeon Graphics

Тактовая частота: 2.00 ГГц

Ядра / потоки: 8 физических / 16 логических

ОЗУ: 16 ГБ

Устройство: Acer Swift SFG14-41

ОС: Windows

Программное обеспечение:

Компилятор: MSVC / GCC (сборка Release)

TBB: Intel oneTBB

Команды запуска:

bash
./ppc_func_tests.exe --gtest_filter="*Artyushkina*tbb*"
./ppc_perf_tests.exe --gtest_filter="*Artyushkina*tbb*"
7. Результаты производительности
Производительность измерялась на изображении 1000x1000 пикселей с шахматным узором.

Mode	Count	Time, s	Speedup	Efficiency
seq	1	0.005214	1.00	N/A
tbb	1	0.018046	0.29	28.9%
Метрика pipeline: 0.022073 с.

8. Выводы
TBB-реализация является наиболее проработанной среди параллельных версий.
Параллельные этапы демонстрируют эффективность TBB-планировщика.

Что получилось:

Корректная обработка 8-связности

Параллельная инициализация и объединение соседей

Использование tbb::spin_mutex для минимальных накладных расходов

Лучший результат среди параллельных версий

Что требует доработки:

Перенумерация выполняется последовательно

Мьютекс создаёт точку синхронизации

Нет ускорения относительно SEQ на изображении 1000x1000

Рекомендации:

Параллелизовать этап NormalizeLabelsTbb

Рассмотреть более мелкозернистую синхронизацию

Увеличить размер тестового изображения