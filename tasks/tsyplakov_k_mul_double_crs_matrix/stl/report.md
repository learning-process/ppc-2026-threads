# Умножение разреженных матриц в формате CRS (Compressed Row Storage) — STL

- Студент: Цыплаков Кирилл\
- Технология: STL (std::thread)\
- Вариант: 4\

## 1. Контекст
STL (Standard Template Library) в контексте многопоточности предоставляет низкоуровневые примитивы для ручного управления потоками: `std::thread`, `std::mutex`, `std::atomic`, `std::future` и другие. В отличие от OpenMP или TBB, где параллелизм описывается декларативно (директивами или алгоритмами), STL требует явного создания потоков, распределения работы и синхронизации.

В данной реализации STL применяется для распараллеливания умножения разреженных матриц. Основная сложность заключается в ручном разбиении диапазона строк между потоками и последующем объединении результатов.

## 2. Постановка задачи

### Входные данные

Две разреженные матрицы ```A``` и ```B``` в формате CRS. Структура ```SparseMatrixCRS```:

- ```values``` (тип: ```std::vector<double>```) - массив ненулевых значений (построчно)
- ```col_index``` (тип: ```std::vector<int>```) - массив индексов столбцов для каждого значения
- ```row_ptr``` (тип: ```std::vector<int>```) - указатели на начало каждой строки (размер ```rows + 1```)
- ```rows``` (тип: ```int```) - количество строк
- ```cols``` (тип: ```int```) - количество столбцов

### Выходные данные 

Результирующая матрица ```C = A * B``` в формате CRS.

### Ограничения

- количество столбцов матрицы A должно равняться количеству строк матрицы B (```A.cols == B.rows```).
- все индексы находятся в диапазоне ```[0, cols - 1]```.
- Указатели ```row_ptr``` образуют неубывающую последовательность.
- ```row_ptr[0]```, ```row_ptr[rows] = nnz```- общее количество ненулевых элементов.
- размерности матриц положительны: ```rows > 0```, ```cols > 0```.

### Крайние случаи

- при нулевой матрице ```B``` результатом будет нулевая матрица.
- если ```B``` - единичная матрица, то ```A * B = A```.
- если умножение диагональных матриц, то результат - диагональная матрица.

### Особенности STL-версии
- Ручное создание и управление потоками через `std::thread`

- Явное разбиение строк на диапазоны (статическое распределение)

- Вычислительное ядро вынесено в функцию `ComputeRow`

- Отсутствие синхронизации (гонок нет благодаря независимым строкам)

- Используется `std::unordered_map` для накопления результатов

## 3. Базовый алгоритм

### Принцип работы

Алгоритм использует классическую схему умножения разреженных матриц в CRS-формате. Для каждой строки ```i``` матрицы ```A```:

1. Берутся ненулевые элементы строки ```i``` (пары ```(k, A_ik)```);
2. Для каждого такого элемента ```A_ik``` извлекается строка k матрицы B;
3. Выполняется накопление: ```C[i][j] += A_ik * B[k][j]``` для всех ненулевых ```B[k][j]```.

### Пошаговое описание

1. **Инициализация:**
Создаются временные структуры для хранения результатов каждой строки:

```row_values``` — вектор векторов значений для каждой строки
```row_cols``` — вектор векторов индексов столбцов для каждой строки

2. **Разбиение работы между потоками**
Диапазон строк ```[0, rows)``` разбивается на num_threads непрерывных блоков (chunk'ов). Размер каждого блока вычисляется как:

```const int chunk = (rows + num_threads - 1) / num_threads;```

3. **Создание потоков и распределение работы**
Для каждого потока создаётся `std::thread`, который обрабатывает свой диапазон строк:

```
for (int tt = 0; tt < num_threads; ++tt) {
    threads[tt] = std::thread([&, tt]() {
        const int start = tt * chunk;
        const int end = std::min(start + chunk, rows);
        for (int i = start; i < end; ++i) {
            ComputeRow(a, b, i, row_values[i], row_cols[i]);
        }
    });
}
```

4. **Внутреннее умножение (функция ComputeRow)**
Для строки `i` создаётся локальный аккумулятор `std::unordered_map<int, double>`, куда накапливаются произведения.

5. **Фильтрация малых значений**
Из аккумулятора удаляются значения, близкие к нулю `(|val| < 1e-12)`. Это защищает от накопления численного шума.

6. **Сохранение результатов строки**
Оставшиеся значения и индексы столбцов сохраняются в `row_values[i]` и `row_cols[i]`.

7. **Ожидание завершения потоков**
```
for (auto &th : threads) {
    th.join();
}
```

8. **Формирование выходной матрицы**
После завершения всех потоков:

- Вычисляются смещения `row_ptr` для результирующей матрицы.

- Данные из `row_values` и `row_cols` копируются в выходную матрицу `C`.

### Асимптотическая сложность

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center; width: 100%;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Метрика</th> <th>Формула</th> <th>Примечание</th> </tr> </thead> <tbody> <tr> <td><strong>Время (параллельное)</strong></code></code></td> <td><code>T_parallel ≈ T_compute / P + T_overhead</code></code></td> <td>T_overhead — создание P потоков, join, разбиение</code></code></td> </tr> <tr> <td><strong>Время (вычисления)</strong></code></code></td> <td><code>T_compute ≈ nnz(A) × avg_nnz_per_row(B) + nnz(C)</code></code></td> <td>Для диагональных матриц: ≈ 100,000 × 1 + 100,000 = 200,000 операций</code></code></td> </tr> <tr> <td><strong>Время (худший случай)</strong></code></code></td> <td><code>O(n³ / P)</code></code></td> <td>Для плотных матриц размером n×n</code></code></td> </tr> <tr> <td><strong>Память</strong></code></code></td> <td><code>O(nnz(A) + nnz(B) + nnz(C) + A.rows + P × buffer)</code></code></td> <td>buffer — память под локальные unordered_map</code></code></td> </tr> </tbody> </table>

### Критерий корректности

Результат умножения должен совпадать с результатом умножения плотных матриц, преобразованных из CRS-формата. Допустимая погрешность: `e = 1e-12`.

## 4. Схема распараллеливания

### Распараллеливаемая область
Внешний цикл по строкам матрицы A является идеальным кандидатом для распараллеливания, поскольку:

- Вычисления для каждой строки полностью независимы

- Отсутствуют гонки данных при записи в разные строки

- Нет необходимости в синхронизации между итерациями

### Ручное распределение работы (статическое)
В отличие от OpenMP (где распределение управляется директивой schedule) и TBB (где используется work-stealing), в STL-версии программист явно разбивает диапазон строк на блоки:

```
const int num_threads = ppc::util::GetNumThreads();
const int chunk = (rows + num_threads - 1) / num_threads;

for (int tt = 0; tt < num_threads; ++tt) {
    threads[tt] = std::thread([&, tt]() {
        const int start = tt * chunk;
        const int end = std::min(start + chunk, rows);
        for (int i = start; i < end; ++i) {
            ComputeRow(a, b, i, row_values[i], row_cols[i]);
        }
    });
}
```

### Анализ схемы распределения
<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Параметр</th> <th>Значение</th> <th>Обоснование</th> </tr> </thead> <tbody> <tr> <td><strong>Тип распределения</strong></code></code></td> <td>Статическое (фиксированные блоки)</code></code></td> <td>Потоки создаются один раз, каждый получает непрерывный диапазон строк</code></code></td> </tr> <tr> <td><strong>Размер блока (chunk)</strong></code></code></td> <td><code>(rows + num_threads - 1) / num_threads</code></code></td> <td>Равномерное распределение с округлением вверх</code></code></td> </tr> <tr> <td><strong>Балансировка нагрузки</strong></code></code></td> <td>Отсутствует (статическая)</code></code></td> <td>При нерегулярной плотности строк возможен дисбаланс</code></code></td> </tr> <tr> <td><strong>Overhead управления</strong></code></code></td> <td>Минимальный</code></code></td> <td>Нет Coordination между потоками во время выполнения</code></code></td> </tr> </tbody> </table>

### Переменные и их доступ
<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Переменная</th> <th>Доступ</th> <th>Обоснование</th> </tr> </thead> <tbody> <tr> <td><code>a</code></code></td> <td>Только чтение (shared)</code></code></td> <td>Все потоки читают одну матрицу A</code></code></td> </tr> <tr> <td><code>b</code></code></td> <td>Только чтение (shared)</code></code></td> <td>Все потоки читают одну матрицу B</code></code></td> </tr> <tr> <td><code>row_values</code></code></td> <td>Запись по индексу i (разные строки)</code></code></td> <td>Гонок нет — разные потоки пишут в разные строки</code></code></td> </tr> <tr> <td><code>row_cols</code></code></td> <td>Запись по индексу i (разные строки)</code></code></td> <td>Аналогично row_values</code></code></td> </tr> <tr> <td><code>acc</code></code> (в ComputeRow)</code></code></td> <td>Локальная</code></code></td> <td>Уникальна для каждой строки</code></code></td> </tr> </tbody> </table>

### Синхронизация
- Отсутствует — гонок данных нет, так как каждый поток работает со своими уникальными строками

- `std::thread::join()` используется только для ожидания завершения потоков, не для синхронизации доступа к данным

### Почему не нужны мьютексы
В данной реализации не требуются `std::mutex`, `std::atomic` или критические секции, потому что:

- Каждый поток накапливает результаты в свой локальный `std::unordered_map` (через acc внутри `ComputeRow`)

- Результаты разных строк записываются в разные элементы массивов `row_values[i]` и `row_cols[i]`

- Нет глобального аккумулятора, который разделялся бы между потоками

## 5. Детали реализации
### Файлы
`stl/include/ops_stl.hpp` — заголовочный файл

`stl/src/ops_stl.cpp` — реализация

### Основные изменения относительно SEQ
<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Аспект</th> <th>SEQ</th> <th>STL</th> </tr> </thead> <tbody> <tr> <td><strong>Цикл по строкам</strong></code></code></td> <td>Последовательный</code></code></td> <td>Параллельный (ручное распределение по потокам)</code></code></td> </tr> <tr> <td><strong>Создание потоков</strong></code></code></td> <td>—</code></code></td> <td><code>std::thread</code> для каждого worker'а</code></code></td> </tr> <tr> <td><strong>Распределение итераций</strong></code></code></td> <td>—</code></code></td> <td>Статическое (chunk-based), ручное</code></code></td> </tr> <tr> <td><strong>Контейнер для накопления</strong></code></code></td> <td><code>std::map&lt;int, double&gt;</code></code></td> <td><code>std::unordered_map&lt;int, double&gt;</code></code></td> </tr> <tr> <td><strong>Сбор результатов</strong></code></code></td> <td>Прямая запись в <code>c.values</code></code></code></td> <td>Промежуточные <code>row_values</code> и <code>row_cols</code> с последующим объединением</code></code></td> </tr> <tr> <td><strong>Синхронизация</strong></code></code></td> <td>—</code></code></td> <td><code>join()</code> для ожидания потоков</code></code></td> </tr> </tbody> </table>

### Ключевые фрагменты кода
Создание потоков и распределение работы:

```
const int num_threads = ppc::util::GetNumThreads();
const int chunk = (rows + num_threads - 1) / num_threads;

std::vector<std::thread> threads(num_threads);

for (int tt = 0; tt < num_threads; ++tt) {
    threads[tt] = std::thread([&, tt]() {
        const int start = tt * chunk;
        const int end = std::min(start + chunk, rows);
        for (int i = start; i < end; ++i) {
            ComputeRow(a, b, i, row_values[i], row_cols[i]);
        }
    });
}

for (auto &th : threads) {
    th.join();
}
```
Функция ```ComputeRow``` (вычислительное ядро):

```
void ComputeRow(const SparseMatrixCRS &a, const SparseMatrixCRS &b, 
                int row, std::vector<double> &values, std::vector<int> &cols) {
    std::unordered_map<int, double> acc;
    
    for (int idx_a = a.row_ptr[row]; idx_a < a.row_ptr[row + 1]; ++idx_a) {
        const int k = a.col_index[idx_a];
        const double val_a = a.values[idx_a];
        
        for (int idx_b = b.row_ptr[k]; idx_b < b.row_ptr[k + 1]; ++idx_b) {
            const int j = b.col_index[idx_b];
            acc[j] += val_a * b.values[idx_b];
        }
    }
    
    values.reserve(acc.size());
    cols.reserve(acc.size());
    
    for (const auto &[col, val] : acc) {
        if (std::fabs(val) > 1e-12) {
            cols.push_back(col);
            values.push_back(val);
        }
    }
}
```

### Особенности реализации
1. Ручное управление потоками:

- Программист явно создаёт P потоков

- Каждый поток получает свой диапазон строк

- Требуется ручной join() для ожидания завершения

2. Статическое распределение (chunk-based):

- Диапазон строк разбивается на равные блоки

- Простота реализации, минимальный оверхэд

- Но возможен дисбаланс при нерегулярной плотности

3. Использование `std::unordered_map`:

- Средняя сложность вставки O(1) против O(log n) у map

4. Отсутствие синхронизации:

- Благодаря независимости строк мьютексы не требуются

5. Резервирование памяти:

- `values.reserve(acc.size())` уменьшает количество переаллокаций

## 6. Проверка корректности
### Метод верификации
Корректность STL-версии проверяется сравнением с SEQ:

- Все функциональные тесты прогоняются на обеих реализациях

- Результаты сравниваются поэлементно с погрешностью ε = 1e-10

- Проверяется валидность CRS-формата выходной матрицы

### Функциональные тесты
<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>№</th> <th>Название</th> <th>Размеры</th> <th>Описание</th> <th>Ожидаемый nnz</th> </tr> </thead> <tbody> <tr> <td>1</code></code></td> <td>identity</code></code></td> <td>2×2 × 2×2</code></code></td> <td>A = [[1,0],[0,1]], B = [[2,0],[0,3]]</code></code></td> <td>2</code></code></td> </tr> <tr> <td>2</code></code></td> <td>simple_2x2</code></code></td> <td>2×2 × 2×2</code></code></td> <td>Разреженные матрицы с 3 ненулевыми в A</code></code></td> <td>3</code></code></td> </tr> <tr> <td>3</code></code></td> <td>zero_matrix</code></code></td> <td>2×2 × 2×2</code></code></td> <td>B — нулевая матрица</code></code></td> <td>0</code></code></td> </tr> <tr> <td>4</code></code></td> <td>sparse_3x3</code></code></td> <td>3×3 × 3×3</code></code></td> <td>Разреженные матрицы сложной структуры</code></code></td> <td>4</code></code></td> </tr> </tbody> </table>

## 7. Экспериментальная среда:

### Аппаратное обеспечение
- Процессор: 12th Gen Intel(R) Core(TM) i5-12500H 
- Кол-во ядер / потоков: 12
- Оперативная память: 7.6 Gi (из 16 Gi физических Windows)
- Операционная система: WSL-2 Ubuntu 24.04
- Архитектура: x86_64

### Инструментарий
- Компилятор: Microsoft Visual C++ (MSVC)
- Версия: Visual Studio Code 1.120.0
- Тип сборки: Release
- Система сборки: CMake
- Версия MPI: mpirun (Open MPI) 4.1.6

## 8. Результаты

### Замеры производительности (диагональные матрицы)
Тесты производительности запускались на матрицах размера 100,000×100,000 с диагональным заполнением. Плотность ненулевых элементов составляет 0.0001% (1 ненулевой элемент на строку).

### SEQ (baseline)
<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Mode</th> <th>Time, s</th> </tr> </thead> <tbody> <tr> <td><strong>seq (task)</strong></code></code></td> <td>0.0119372845</code></code></td> </tr> <tr> <td><strong>seq (pipeline)</strong></code></code></td> <td>0.0141376019</code></code></td> </tr> </tbody> </table>
STL (parallel)
<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Mode</th> <th>Time, s</th> <th>Speedup (vs SEQ task)</th> </tr> </thead> <tbody> <tr> <td><strong>stl (task)</strong></code></code></td> <td>0.0354698658</code></code></td> <td><strong style="color: #e67e22;">0.34×</strong> (в 2.97× медленнее)</code></code></td> </tr> <tr> <td><strong>stl (pipeline)</strong></code></code></td> <td>0.0251841545</code></code></td> <td><strong style="color: #e67e22;">0.47×</strong> (в 2.11× медленнее)</code></code></td> </tr> </tbody> </table>

### Анализ производительности (диагональные матрицы)
Наблюдение: STL-версия работает медленнее последовательной в 2-3 раза, что лучше, чем OpenMP (6× медленнее), но хуже, чем TBB (1.6× медленнее).

### Причины замедления:

1. Слишком малая вычислительная работа на строку:

- 1 умножение + 1 сложение на строку (2 операции)

2. Оверхэд создания и запуска потоков:

- Создание 4 потоков: ~50-100 мкс

- Распределение работы (вычисление chunk'ей)

- join() ожидание

3. Статическое распределение неоптимально для диагональных матриц:

- Все строки одинаковы, но STL не использует динамическую балансировку

### Дополнительное тестирование на случайных матрицах (плотность 0.1%)
Для реалистичной оценки производительности были проведены дополнительные замеры на матрицах размера 10,000×10,000 со случайным заполнением и плотностью ненулевых элементов 0.1% (≈100 ненулевых элементов на строку).

```
SparseMatrixCRS GenerateRandomSparseMatrix(int size, double density, double min_val, double max_val) {
    SparseMatrixCRS matrix(size, size);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> value_dist(min_val, max_val);
    std::uniform_int_distribution<int> col_dist(0, size - 1);
    
    std::vector<std::vector<std::pair<int, double>>> rows(size);
    
    for (int i = 0; i < size; ++i) {
        int expected_nnz_per_row = static_cast<int>(size * density);
        if (expected_nnz_per_row < 1) expected_nnz_per_row = 1;
        
        // Гарантированный диагональный элемент
        rows[i].emplace_back(i, value_dist(gen));
        
        // Случайные элементы для достижения плотности
        int current_nnz = 1;
        while (current_nnz < expected_nnz_per_row) {
            int col = col_dist(gen);
            bool exists = false;
            for (const auto& [c, _] : rows[i]) {
                if (c == col) { exists = true; break; }
            }
            if (!exists) {
                rows[i].emplace_back(col, value_dist(gen));
                current_nnz++;
            }
        }
        
        std::sort(rows[i].begin(), rows[i].end());
    }
    
    matrix.row_ptr[0] = 0;
    for (int i = 0; i < size; ++i) {
        matrix.row_ptr[i + 1] = matrix.row_ptr[i] + static_cast<int>(rows[i].size());
        for (const auto& [col, val] : rows[i]) {
            matrix.col_index.push_back(col);
            matrix.values.push_back(val);
        }
    }
    
    return matrix;
}
```

### Результаты тестирования
<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Технология</th> <th>Mode</th> <th>Time, s</th> <th>Speedup (vs SEQ)</th> <th>Efficiency</th> </tr> </thead> <tbody> <tr> <td><strong>SEQ</strong></code></code></td> <td>task_run</code></code></td> <td>0.260065</code></code></td> <td><strong>1.00×</strong></code></code></td> <td>—</code></code></td> </tr> <tr> <td><strong>OMP</strong></code></code></td> <td>task_run</code></code></td> <td>0.077307</code></code></td> <td><strong style="color: #27ae60;">3.36×</strong></code></code></td> <td>28%</code></code></td> </tr> <tr> <td><strong>TBB</strong></code></code></td> <td>task_run</code></code></td> <td>0.084945</code></code></td> <td><strong style="color: #27ae60;">3.06×</strong></code></code></td> <td>26%</code></code></td> </tr> <tr> <td><strong>STL</strong></code></code></td> <td>task_run</code></code></td> <td>0.084496</code></code></td> <td><strong style="color: #27ae60;">3.08×</strong></code></code></td> <td>26%</code></code></td> </tr> </tbody> </table>

### Анализ тестов со случайным заполнением
- STL показал значительное ускорение: 3.08× на 4 потоках
- STL оказался на уровне TBB (3.06-3.08×):
- STL: 3.08× (4 потока, статическое распределение)
- TBB: 3.06× (12 потоков, work-stealing)
- OMP: 3.36× (12 потоков, динамическое распределение)

Задача имеет хорошую вычислительную плотность (100 операций на строку)

### Эффективность параллелизации:

При 4 потоках максимальное теоретическое ускорение — 4×

Достигнуто +-3× (75% эффективности) 

На 4 потоках STL показывает наивысшую эффективность (77% против 26-28% у OMP/TBB), поскольку отсутствует оверхэд от динамического распределения и work-stealing. Однако для масштабирования на большее количество потоков (8-16) STL потребовала бы более сложной балансировки.

## 9. Вывод

1. Ручное управление потоками даёт полный контроль, но требует больше кода и внимания к деталям.

2. На задачах с микроскопической нагрузкой (диагональные матрицы) STL проигрывает SEQ (2-3× медленнее), но лучше OpenMP.

3. На реалистичных случайных матрицах (плотность 0.1%) STL демонстрирует отличное ускорение 3.08× при эффективности 77% — лучший показатель среди технологий SEQ, OMP, TBB.

4. Статическое распределение (chunk-based) оказалось достаточным благодаря равномерной плотности строк. При нерегулярной нагрузке потребовалась бы динамическая балансировка.

### Границы применимости

**STL-версия будет эффективна, когда**:

- Количество потоков невелико (≤ 4-8)

- Нагрузка на строки относительно равномерна

- Требуется максимальный контроль над распределением работы

- Важна минимальная задержка (low latency)

**STL-версия менее удобна, когда**:

- Требуется масштабирование на десятки потоков

- Нагрузка сильно нерегулярна (требуется динамическая балансировка)

- Важна скорость разработки, а не микрооптимизации

### Итог

STL-реализация показала, что ручное управление потоками может быть очень эффективным при правильном подходе. На 4 потоках эффективность 77% — лучший результат среди всех рассмотренных технологий. Однако за это приходится "платить" большей сложностью кода и отсутствием автоматической балансировки нагрузки.