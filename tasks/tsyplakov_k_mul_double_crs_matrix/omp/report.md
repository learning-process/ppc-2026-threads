# Умножение разреженных матриц в формате CRS (Compressed Row Storage) — OMP

- Студент: Цыплаков Кирилл
- Технология: OMP
- Вариант: 4

## 1. Контекст

OpenMP (Open Multi-Processing) — стандарт параллельного программирования для систем с общей памятью, основанный на использовании директив компилятора. В данной реализации OpenMP применяется для распараллеливания внешнего цикла по строкам матрицы A, поскольку вычисления для разных строк независимы.

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

### Особенности OMP-версии

- Распараллеливается обработка строк результирующей матрицы C
- Каждый поток обрабатывает свой набор строк независимо
- Используется ```std::unordered_map``` для накопления результатов (вместо ```std::map``` в SEQ)

## 3. Базовый алгоритм

### Принцип работы

Алгоритм использует классическую схему умножения разреженных матриц в CRS-формате. Для каждой строки `i` матрицы `A`:

- Берутся ненулевые элементы строки `i` (пары `(k, A_ik)`)
- Для каждого такого элемента `A_ik` извлекается строка `k` матрицы `B`
- Выполняется накопление: `C[i][j] += A_ik * B[k][j]` для всех ненулевых `B[k][j]`

### Математическая основа

Умножение матриц `C = A × B` определяется формулой:

`C[i][j] = Σ_{k=0}^{n-1} A[i][k] × B[k][j]`

В разреженном формате CRS итерация идёт только по ненулевым элементам, что значительно ускоряет вычисления при малой плотности заполнения.

### Пошаговое описание

#### Шаг 1: Инициализация

Создаются временные структуры для хранения результатов каждой строки:

`row_values` — вектор векторов значений для каждой строки

row_cols — вектор векторов индексов столбцов для каждой строки

#### Шаг 2: Параллельный обход строк матрицы A

Для каждой строки `i` от `0` до `A.rows-1` выполняется:

`std::unordered_map<int, double> acc;`

#### Шаг 3: Внутреннее умножение

Для каждого ненулевого элемента `(k, A_ik)` из строки `i` матрицы `A`:

Из матрицы `B` извлекается строка `k` (все ненулевые элементы в этой строке)

Для каждого элемента `(j, B_kj)` из строки k:

`acc[j] += A_ik × B_kj`

#### Шаг 4: Фильтрация малых значений

Из `acc` удаляются значения, близкие к нулю `(|val| < 1e-12)`. Это защищает от накопления численного шума.

#### Шаг 5: Сохранение результатов строки

Оставшиеся значения и индексы столбцов сохраняются:
`
row_cols[i].push_back(col);
row_values[i].push_back(val);
`

#### Шаг 6: Формирование выходной матрицы

После завершения всех параллельных итераций:

- Вычисляются смещения `row_ptr` для результирующей матрицы
- Данные из row_values` и `row_cols` копируются в выходную матрицу C

### Асимптотическая сложность

Обозначения:

`nnz(A)` — количество ненулевых элементов в матрице A
`nnz(B)` — количество ненулевых элементов в матрице B
`nnz(C)` — количество ненулевых элементов в результирующей матрице C

`n` — размерность матриц (для квадратного случая)
`P` — количество потоков OpenMP

### Время выполнения

`T_parallel = T_compute / P + T_overhead`, где:

`T_compute = nnz(A) × avg_nnz_per_row(B) + nnz(C)`
`T_overhead` — затраты на создание потоков, распределение итераций, синхронизацию

В худшем случае (плотные матрицы): `O(n^3 / P)`
В разреженном случае (текущая задача): `O(100000 / P + overhead)`

Память:

`O(nnz(A) + nnz(B) + nnz(C) + A.rows + P × buffer_size)`
где `buffer_size` — память под локальные unordered_map в каждом потоке

### Критерий корректности

Результат умножения должен совпадать с результатом умножения плотных матриц, преобразованных из CRS-формата. Допустимая погрешность: `e = 1e-12` для сравнения значений с плавающей точкой.

## 4. Схема распараллеливания

### Распараллеливаемая область

Внешний цикл по строкам матрицы `A` `(for (int i = 0; i < rows; ++i))` является идеальным кандидатом для распараллеливания, поскольку:

- Вычисления для каждой строки `i` полностью независимы
- Отсутствуют гонки данных при записи в разные строки результирующей матрицы
- Нет необходимости в синхронизации между итерациями

### Директивы OpenMP

`#pragma omp parallel for schedule(dynamic) default(none) shared(a, b, row_values, row_cols, rows)
for (int i = 0; i < rows; ++i) {
    // обработка строки i
}
`

### Анализ параметров директивы

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; font-family: Arial, sans-serif; width: 100%;">
  <caption style="caption-side: top; margin-bottom: 8px; font-weight: bold; text-align: center;">Анализ параметров директивы OpenMP</caption>
  <thead>
    <tr style="background-color: #f2f2f2;">
      <th>Параметр</th>
      <th>Значение</th>
      <th>Обоснование</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>parallel for</code></code></td>
      <td>Распараллеливание цикла</code></code></td>
      <td>Создаётся команда потоков, итерации распределяются между ними</code></code></td>
    </tr>
    <tr>
      <td><code>schedule(dynamic)</code></code></td>
      <td>Динамическое распределение итераций</code></code></td>
      <td>Строки могут иметь разное количество ненулевых элементов; dynamic позволяет балансировать нагрузку</code></code></td>
    </tr>
    <tr>
      <td><code>default(none)</code></code></td>
      <td>Явное указание атрибутов переменных</code></code></td>
      <td>Требует перечислить все переменные с атрибутами shared/private, что повышает безопасность</code></code></td>
    </tr>
    <tr>
      <td><code>shared(a, b, row_values, row_cols, rows)</code></code></td>
      <td>Общие переменные</code></code></td>
      <td>Все потоки читают одни и те же матрицы A и B; row_values/row_cols — массивы для сбора результатов</code></code></td>
    </tr>
  </tbody>
</table>

### Переменные и их атрибуты

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; font-family: Arial, sans-serif; width: 100%;">
  <caption style="caption-side: top; margin-bottom: 8px; font-weight: bold; text-align: center;">Переменные и их атрибуты</caption>
  <thead>
    <tr style="background-color: #f2f2f2;">
      <th>Переменная</th>
      <th>Атрибут</th>
      <th>Обоснование</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>a</code></code></td>
      <td><code>shared</code></code></td>
      <td>Только чтение, все потоки обращаются к одной матрице</code></code></td>
    </tr>
    <tr>
      <td><code>b</code></code></td>
      <td><code>shared</code></code></td>
      <td>Только чтение, все потоки обращаются к одной матрице</code></code></td>
    </tr>
    <tr>
      <td><code>rows</code></code></td>
      <td><code>shared</code></code></td>
      <td>Константа, только чтение</code></code></td>
    </tr>
    <tr>
      <td><code>row_values</code></code></td>
      <td><code>shared</code></code></td>
      <td>Доступ по индексу i — разные потоки работают с разными строками</code></code></td>
    </tr>
    <tr>
      <td><code>row_cols</code></code></td>
      <td><code>shared</code></code></td>
      <td>Аналогично row_values</code></code></td>
    </tr>
    <tr>
      <td><code>i</code></code></td>
      <td><code>private</code> (неявно)</code></code></td>
      <td>Каждый поток имеет свою копию итератора</code></code></td>
    </tr>
    <tr>
      <td><code>acc</code></code></td>
      <td><code>private</code> (локальная)</code></code></td>
      <td>Создаётся внутри цикла, уникальна для каждой итерации</code></code></td>
    </tr>
  </tbody>
</table>

### Отсутствие редукции и критических секций

В данной реализации не требуются ```reduction```, ```atomic``` или `critical`, потому что:

- Каждый поток накапливает результаты в свой локальный `std::unordered_map (acc)`
- Результаты разных строк записываются в разные элементы массивов `row_values[i]` и `row_cols[i]`
- Нет глобального аккумулятора, который разделялся бы между потоками

### Барьеры и синхронизация

В конце `parallel for` существует неявный барьер — все потоки ожидают завершения обработки своих строк

Дополнительных явных барьеров не требуется.

## 5. Детали реализации

### Файлы

`omp/include/ops_omp.hpp` — заголовочный файл\
`omp/src/ops_omp.cpp` — реализация

### Основные изменения относительно SEQ

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; font-family: Arial, sans-serif; width: 100%;">
  <caption style="caption-side: top; margin-bottom: 8px; font-weight: bold; text-align: center;">Основные изменения относительно SEQ</caption>
  <thead>
    <tr style="background-color: #f2f2f2;">
      <th>Аспект</th>
      <th>SEQ</th>
      <th>OMP</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>Цикл по строкам</strong></td>
      <td>Последовательный</td>
      <td>Параллельный (<code>#pragma omp parallel for</code>)</td>
    </tr>
    <tr>
      <td><strong>Контейнер для накопления</strong></td>
      <td><code>std::map&lt;int, double&gt;</code></td>
      <td><code>std::unordered_map&lt;int, double&gt;</code></td>
    </tr>
    <tr>
      <td><strong>Сбор результатов</strong></td>
      <td>Прямая запись в <code>c.values</code></td>
      <td>Промежуточные <code>row_values</code> и <code>row_cols</code> с последующим объединением</td>
    </tr>
    <tr>
      <td><strong>Распределение итераций</strong></td>
      <td>—</td>
      <td><code>schedule(dynamic)</code></td>
    </tr>
  </tbody>
</table>

### Ключевые фрагменты кода

**Параллельный цикл:**

```
#pragma omp parallel for schedule(dynamic) default(none) shared(a, b, row_values, row_cols, rows)
for (int i = 0; i < rows; ++i) {
    std::unordered_map<int, double> acc;  // локальный для каждого потока
    
    for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
        int k = a.col_index[j];
        double val_a = a.values[j];
        
        for (int zz = b.row_ptr[k]; zz < b.row_ptr[k + 1]; ++zz) {
            int j1 = b.col_index[zz];
            acc[j1] += val_a * b.values[zz];  // локальная запись, без гонок
        }
    }
    
    // Сохранение результатов для строки i
    row_values[i].reserve(acc.size());
    row_cols[i].reserve(acc.size());
    
    for (const auto &[col, val] : acc) {
        if (std::fabs(val) > 1e-12) {
            row_cols[i].push_back(col);
            row_values[i].push_back(val);
        }
    }
}
```

Сбор результатов в финальную матрицу:

```
SparseMatrixCRS c(a.rows, b.cols);

// Вычисление смещений row_ptr
for (int i = 0; i < c.rows; ++i) {
    c.row_ptr[i + 1] = c.row_ptr[i] + static_cast<int>(row_values[i].size());
}

// Копирование данных
for (int i = 0; i < c.rows; ++i) {
    c.values.insert(c.values.end(), row_values[i].begin(), row_values[i].end());
    c.col_index.insert(c.col_index.end(), row_cols[i].begin(), row_cols[i].end());
}
```

### Особенности реализации

- Использование `std::unordered_map` вместо `std::map`:
- `unordered_map` имеет среднюю сложность вставки `O(1)` против `O(log n)` у `map`. Это уменьшает оверхэд в параллельной версии

### Промежуточное хранение результатов

1. Каждая строка сначала накапливается в локальный ```unordered_map```
2. Затем результаты сохраняются в `row_values[i]` и `row_cols[i]`

Финальная сборка выполняется последовательно вне параллельной области

### Резервирование памяти

`row_values[i].reserve(acc.size())` — уменьшает количество переаллокаций

### Schedule(dynamic)

Выбран потому, что строки могут иметь разную «плотность» ненулевых элементов
Динамическое распределение позволяет сбалансировать нагрузку между потоками

## 6. Проверка корректности

### Метод верификации

Корректность OMP-версии проверяется сравнением с SEQ:

- Все функциональные тесты (4 тестовых случая) прогоняются на обеих реализациях
- Результаты сравниваются поэлементно с допустимой погрешностью ε = 1e-10
- Дополнительно проверяется валидность CRS-формата выходной матрицы

### Функциональные тесты

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>№</th> <th>Название</th> <th>Размеры</th> <th>Описание</th> <th>Ожидаемый nnz</th> </tr> </thead> <tbody> <tr> <td>1</td> <td>identity</td> <td>2×2 × 2×2</td> <td>A = [[1,0],[0,1]], B = [[2,0],[0,3]]</td> <td>2</td> </tr> <tr> <td>2</td> <td>simple_2x2</td> <td>2×2 × 2×2</td> <td>Разреженные матрицы с 3 ненулевыми в A</td> <td>3</td> </tr> <tr> <td>3</td> <td>zero_matrix</td> <td>2×2 × 2×2</td> <td>B — нулевая матрица</td> <td>0</td> </tr> <tr> <td>4</td> <td>sparse_3x3</td> <td>3×3 × 3×3</td> <td>Разреженные матрицы сложной структуры</td> <td>4</td> </tr> </tbody> </table>

## 7. Экспериментальная среда

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

### Замеры производительности

Тесты производительности запускались на матрицах размера `100000×100000` с диагональным заполнением. Плотность ненулевых элементов составляет 0.0001%.

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Mode</th> <th>Workers</th> <th>Time, s</th> <th>Speedup (vs SEQ)</th> <th>Efficiency</th> </tr> </thead> <tbody> <tr> <td><strong>omp (task)</strong></td> <td>12 (по умолчанию)</td> <td>0.0717499506</td> <td>0.162</td> <td>1.35%</td> </tr> <tr> <td><strong>omp (pipeline)</strong></td> <td>12 (по умолчанию)</td> <td>0.0660709260</td> <td>0.176</td> <td>1.47%</td> </tr> </tbody> </table>

### Сравнение с SEQ

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;">
  <thead>
    <tr style="background-color: #f2f2f2;">
      <th>Показатель</th>
      <th>SEQ (task)</th>
      <th>OMP (task)</th>
      <th>OMP (pipeline)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Время, с</code></code></td>
      <td>0.0116648674</code></code></td>
      <td>0.0717499506</code></code></td>
      <td>0.0660709260</code></code></td>
    </tr>
    <tr>
      <td>Относительное время</code></code></td>
      <td>1×</code></code></td>
      <td>6.15× медленнее</code></code></td>
      <td>5.66× медленнее</code></code></td>
    </tr>
  </tbody>
</table>

**Причины замедления**:

1. Слишком малая вычислительная работа на строку:

- Матрицы диагональные, поэтому каждая строка A содержит ровно 1 ненулевой элемент
- Каждая строка B также содержит 1 ненулевой элемент
- На строку выполняется всего 1 умножение + 1 сложение

2. Оверхэд OpenMP превышает выигрыш:

- Создание 12 потоков: 50-100 мкс
- Динамическое распределение итераций: дополнительные накладные расходы
- Синхронизация (неявный барьер): ещё 10-50 мкс
- На 100000 итераций суммарный оверхэд становится доминирующим

3. Промежуточное хранение результатов:

- Векторы `row_values` и `row_cols` требуют дополнительного копирования
- В SEQ данные записываются напрямую в выходную матрицу

4. ```schedule(dynamic)``` неоптимален:

- Для диагональных матриц все строки одинаковы (1 ненулевой элемент)
- `schedule(static)` был бы эффективнее, так как не требует coordination

### Дополнительное тестирование на случайных матрицах

Для того, чтобы получить основное преимущество технологии OMP - ускорение, матрицы должны быть менее разреженными. Для этого в тесты добавим случайное заполнение матриц:

```
...хэдеры...

namespace tsyplakov_k_mul_double_crs_matrix {

class TsyplakovKRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kSize_ = 10000;           // матрица размером 10000 на 10000
  const double kDensity_ = 0.001;      // 0.1% плотности (100 ненулевых на строку)
  const double kMinValue_ = -10.0;     // минимальное значение ненулевого элемента
  const double kMaxValue_ = 10.0;      // максимальное значение ненулевого элемента

  InType input_data_{};

  // Генерация случайной разреженной матрицы в формате CRS
  SparseMatrixCRS GenerateRandomSparseMatrix(int size, double density, double min_val, double max_val) {
    SparseMatrixCRS matrix(size, size);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> density_dist(0.0, 1.0);
    std::uniform_real_distribution<double> value_dist(min_val, max_val);
    std::uniform_int_distribution<int> col_dist(0, size - 1);
    
    std::vector<std::vector<std::pair<int, double>>> rows(size);
    
    for (int i = 0; i < size; ++i) {
      // Ожидаемое количество ненулевых элементов в строке i
      int expected_nnz_per_row = static_cast<int>(size * density);
      if (expected_nnz_per_row < 1) expected_nnz_per_row = 1;
      
      // Используем оба подхода для более реалистичного распределения:
      // 1) Часть элементов — случайные позиции
      // 2) Часть элементов — гарантированные (чтобы не было пустых строк)
      
      // Добавляем хотя бы 1 гарантированный элемент на строку (диагональный)
      rows[i].emplace_back(i, value_dist(gen));
      
      // Добавляем случайные элементы для достижения желаемой плотности
      int current_nnz = 1;
      while (current_nnz < expected_nnz_per_row) {
        int col = col_dist(gen);
        // Избегаем дублирования столбцов в одной строке
        bool exists = false;
        for (const auto& [c, _] : rows[i]) {
          if (c == col) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          rows[i].emplace_back(col, value_dist(gen));
          current_nnz++;
        }
      }
      
      // Сортируем по индексу столбца (требование CRS формата)
      std::sort(rows[i].begin(), rows[i].end());
    }
    
    // Заполняем структуру CRS
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

  void SetUp() override {
    // Генерируем две случайные разреженные матрицы с одинаковой плотностью
    SparseMatrixCRS a = GenerateRandomSparseMatrix(kSize_, kDensity_, kMinValue_, kMaxValue_);
    SparseMatrixCRS b = GenerateRandomSparseMatrix(kSize_, kDensity_, kMinValue_, kMaxValue_);
    
    input_data_ = {.a = a, .b = b};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Проверяем, что размерность выходной матрицы корректна
    if (output_data.rows != kSize_ || output_data.cols != kSize_) {
      return false;
    }
    
    // Проверяем, что row_ptr имеет правильный размер
    if (static_cast<int>(output_data.row_ptr.size()) != kSize_ + 1) {
      return false;
    }
    
    // Проверяем, что row_ptr является неубывающей последовательностью
    for (int i = 0; i < kSize_; ++i) {
      if (output_data.row_ptr[i] > output_data.row_ptr[i + 1]) {
        return false;
      }
    }
    
    // Проверяем, что values и col_index имеют одинаковый размер
    if (output_data.values.size() != output_data.col_index.size()) {
      return false;
    }
    
    // Пустой результат допустим, если произведение дало нулевую матрицу
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TsyplakovKRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TsyplakovKTestTaskOMP>(PPC_SETTINGS_tsyplakov_k_mul_double_crs_matrix);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = TsyplakovKRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TsyplakovKRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tsyplakov_k_mul_double_crs_matrix
```

После проведения тестирования были получены следующие результаты:

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Технология</th> <th>Mode</th> <th>Time, s</th> <th>Speedup (vs SEQ)</th> <th>Efficiency</th> </tr> </thead> <tbody> <tr> <td><strong>SEQ</strong></code></code></td> <td>task_run</code></code></td> <td>0.260065</code></code></td> <td><strong>1.00×</strong></code></code></td> <td>—</code></code></td> </tr> <tr> <td><strong>OMP</strong></code></code></td> <td>task_run</code></code></td> <td>0.077307</code></code></td> <td><strong style="color: #27ae60;">3.36×</strong></code></code></td> <td>28%</code></code></td> </tr> <tr> <td><strong>TBB</strong></code></code>（для сравнения）</code></code></td> <td>task_run</code></code></td> <td>0.084945</code></code></td> <td><strong style="color: #27ae60;">3.06×</strong></code></code></td> <td>26%</code></code></td> </tr> </tbody> <table>

### Анализ тестирования с матрицами со случайной генерацией

1. OpenMP показал значительное ускорение: 3.36× на 12 потоках\
2. OpenMP оказался немного быстрее TBB (на ~10%):

- Меньший оверхэд на директивах компилятора

- Отсутствие промежуточного хранения результатов

- Более прямолинейная работа с памятью

3. Эффективность параллелизации:

При 12 потоках максимальное теоретическое ускорение — 12×

Достигнуто ~3× из-за последовательной части алгоритма (сбор результатов, накладные расходы)

## 9. Выводы

OpenMP — мощный инструмент для задач с «тяжёлыми» итерациями, где вычислительная работа на итерацию значительно превышает оверхэд.

На задачах с микроскопической нагрузкой на итерацию (как в данном случае — 1 операция) OpenMP неэффективен и может приводить к замедлению.

Выбор стратегии распараллеливания критически важен:

```schedule(dynamic)``` хорош для несбалансированных нагрузок
`schedule(static)` лучше для однородных задач

Для очень лёгких итераций стоит использовать schedule(static, chunk_size) с блочным разбиением

### Границы применимости

OpenMP-версия будет эффективна, когда:

- На одну итерацию приходится ≥ 100-1000 арифметических операций
- Матрицы имеют плотность >= 0.1% или размер > 1,000,000
- Используется статическое распределение для однородных задач
- На текущей конфигурации (диагональные матрицы 100000×100000) SEQ остаётся лучшим выбором.
