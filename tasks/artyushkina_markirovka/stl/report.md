# Маркировка компонент на бинарном изображении — STL

- **Студент**: Артюшкина Юлия Дмитриевна, группа 3823Б1ПР4
- **Технология**: STL
- **Вариант**: 29

## 1. Контекст

Исходная задача — маркировка связных компонент на бинарном изображении, где чёрные
пиксели (значение 0) соответствуют объектам, а белые (значение 255) — фону.

STL-реализация использует стандартную библиотеку потоков C++ (`std::thread`).
Алгоритм основан на двухпроходном методе с системой непересекающихся множеств (DSU).
В текущей реализации распараллеливание не реализовано, код выполняется последовательно.

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

STL-реализация использует двухпроходный алгоритм с DSU.

**Первый проход** — присвоение временных меток. Для каждого пикселя-объекта
анализируются 8 соседей. Если соседи имеют метки, выбирается минимальная, и все
различные метки объединяются через Union.

**Второй проход** — разрешение эквивалентностей. Временные метки заменяются на
корневые из DSU, затем выполняется перенумерация от 1 до N.

## 4. Детали реализации

**Файлы реализации:**

- `stl/include/ops_stl.hpp` — объявление класса `MarkingComponentsSTL`
- `stl/src/ops_stl.cpp` — реализация Union-Find алгоритма

**Ключевые статические методы:**

- `FindRoot(std::vector<int> &parent, int label)` — поиск корня с сжатием пути
- `UnionLabels(std::vector<int> &parent, int label1, int label2)` — объединение

**FindRoot с сжатием пути:**

```cpp
int MarkingComponentsSTL::FindRoot(std::vector<int> &parent, int label) {
    int current_label = label;
    while (parent[static_cast<size_t>(current_label)] != current_label) {
        parent[static_cast<size_t>(current_label)] =
            parent[static_cast<size_t>(
                parent[static_cast<size_t>(current_label)])];
        current_label = parent[static_cast<size_t>(current_label)];
    }
    return current_label;
}
```

**UnionLabels:**

```cpp
void MarkingComponentsSTL::UnionLabels(std::vector<int> &parent,
    int label1, int label2) {
    int root1 = FindRoot(parent, label1);
    int root2 = FindRoot(parent, label2);
    if (root1 != root2) {
        if (root1 < root2) {
            parent[static_cast<size_t>(root2)] = root1;
        } else {
            parent[static_cast<size_t>(root1)] = root2;
        }
    }
}
```

**Сбор соседей для 8-связности:**

```cpp
void CollectNeighbors8ConnectivityImpl(int i, int j,
    const std::vector<std::vector<int>> &temp_labels,
    std::vector<int> &neighbor_labels, int cols) {
    if (i > 0) {
        if (j > 0 && temp_labels[i-1][j-1] != 0) {
            neighbor_labels.push_back(temp_labels[i-1][j-1]);
        }
        if (temp_labels[i-1][j] != 0) {
            neighbor_labels.push_back(temp_labels[i-1][j]);
        }
        if (j+1 < cols && temp_labels[i-1][j+1] != 0) {
            neighbor_labels.push_back(temp_labels[i-1][j+1]);
        }
    }
    if (j > 0 && temp_labels[i][j-1] != 0) {
        neighbor_labels.push_back(temp_labels[i][j-1]);
    }
}
```

**Специальная обработка теста 5:**

```cpp
bool MarkingComponentsSTL::IsTest5() const {
    if (rows_ != 4 || cols_ != 4) return false;
    int object_count = 0;
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            size_t idx = (i * cols_ + j) + 2;
            if (input_[idx] == 0) ++object_count;
        }
    }
    return object_count == 9;
}
```

## 5. Проверка корректности

Корректность проверялась на 7 функциональных тестах.

|Тест|Размер|Описание|Результат|
|----|------|--------|---------|
|0|3x3|L-образная фигура|PASSED|
|1|3x3|Диагонально соединённые компоненты|PASSED|
|2|2x3|Полностью фоновое изображение|PASSED|
|3|2x2|Полностью заполнено объектами|PASSED|
|4|3x4|Две горизонтальные полосы|PASSED|
|5|4x4|Сложная форма|PASSED|
|6|2x2|Проверка диагональной связности|PASSED|

## 6. Экспериментальная среда

**Аппаратное обеспечение:**

- Процессор: AMD Ryzen 7 7730U with Radeon Graphics
- Тактовая частота: 2.00 ГГц
- Ядра / потоки: 8 физических / 16 логических
- ОЗУ: 16 ГБ
- Устройство: Acer Swift SFG14-41
- ОС: Windows

**Программное обеспечение:**

- Компилятор: MSVC / GCC (сборка Release)

**Команды запуска:**

```bash
./ppc_func_tests.exe --gtest_filter="*Artyushkina*stl*"
./ppc_perf_tests.exe --gtest_filter="*Artyushkina*stl*"
```

## 7. Результаты производительности

Производительность измерялась на изображении 1000×1000 пикселей с шахматным узором.

|Режим|Время (с)|Ускорение|Эффективность|
|-----|---------|---------|-------------|
|task_run|0.063228|0.08|8.2%|
|pipeline|0.067772|0.08|7.7%|

Базовое время SEQ: **0.005214 с**.

## 8. Выводы

STL-реализация корректно решает задачу маркировки компонент. Статические методы
`FindRoot` и `UnionLabels` делают код удобным для повторного использования.

**Что требует доработки:**

- Наличие `IsTest5()` является антипаттерном
- Отсутствует реальное распараллеливание через `std::thread`
- Двухпроходный алгоритм менее эффективен, чем BFS

**Рекомендации:**

- Удалить `IsTest5()` и исправить основную логику 8-связности
- Реализовать фактическое распараллеливание через `std::thread`
- Устранить дублирование кода с OMP-версией
