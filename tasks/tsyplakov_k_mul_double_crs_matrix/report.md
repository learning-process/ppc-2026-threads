# Умножение разреженных матриц в формате CRS (Compressed Row Storage)

- Студент: Цыплаков Кирилл, группа 3823Б1ПР2

- Вариант: 4

- Local reports: seq/report.md, omp/report.md, tbb/report.md, stl/report.md, all/report.md

## 1. Введение
Умножение разреженных матриц является одной из ключевых операций в вычислительной линейной алгебре. Оно широко применяется в методах конечных элементов, решении систем линейных уравнений, обработке графов, машинном обучении и других областях, где возникают матрицы большой размерности с преимущественно нулевыми элементами.

Данная задача особенно подходит для сравнения разных моделей параллелизма по следующим причинам:

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%; font-family: Arial, sans-serif;">
  <thead>
    <tr style="background-color: #f2f2f2;">
      <th>Характеристика</th>
      <th>Описание</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>Высокая степень независимости</strong></td>
      <td>Вычисления для разных строк результирующей матрицы полностью независимы</code></code></td>
    </tr>
    <tr>
      <td><strong>Неравномерная нагрузка</strong></td>
      <td>Разреженная структура создаёт неравномерное распределение вычислительной работы</code></code></td>
    </tr>
    <tr>
      <td><strong>Различные стратегии балансировки</strong></td>
      <td>Можно применить статическое и динамическое распределение итераций</code></code></td>
    </tr>
    <tr>
      <td><strong>Масштабируемость</strong></td>
      <td>Задача хорошо параллелится как на уровне одного узла, так и на кластере</code></code></td>
    </tr>
  </tbody>
</table>

В работе реализованы и сравнены пять вариантов:

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%; font-family: Arial, sans-serif;">
  <thead>
    <tr style="background-color: #f2f2f2;">
      <th>Технология</th>
      <th>Описание</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>SEQ</strong></code></code></code></code></td>
      <td>Последовательная эталонная версия</code></code></code></code></td>
    </tr>
    <tr>
      <td><strong>OMP</strong></code></code></code></code></td>
      <td>Распараллеливание с помощью OpenMP (директивы компилятора)</code></code></code></code></td>
    </tr>
    <tr>
      <td><strong>TBB</strong></code></code></code></code></td>
      <td>Распараллеливание с помощью Intel oneTBB (задачи и work-stealing)</code></code></code></code></td>
    </tr>
    </tr>
      <td><strong>STL</strong></code></code></code></code></td>
      <td>Ручное управление потоками через <code>std::thread</code></code></code></code></code></td>
    </tr>
    <tr>
      <td><strong>ALL</strong></code></code></code></code></td>
      <td>Гибридная версия (MPI + OpenMP) — двухуровневый параллелизм</code></code></code></code></td>
    </tr>
  </tbody>
</table>

## 2. Единая постановка задачи

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

## 3. Единая методика эксперимента

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

### Переменные окружения

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%; font-family: Arial, sans-serif;">
  <thead>
    <tr style="background-color: #f2f2f2;">
      <th>Технология</th>
      <th>Переменные окружения</th>
      <th>Пример</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>SEQ</strong></code></code></td>
      <td>(не требуются)</code></code></td>
      <td>—</code></code></td>
    </tr>
    <tr>
      <td><strong>OMP</strong></code></code></td>
      <td><code>OMP_NUM_THREADS=N</code></code></code></code></td>
      <td><code>export OMP_NUM_THREADS=12</code></code></code></code></td>
    </tr>
    <tr>
      <td><strong>TBB</strong></code></code></td>
      <td><code>PPC_NUM_THREADS=N</code></code></code></code></td>
      <td><code>export PPC_NUM_THREADS=12</code></code></code></code></td>
    </tr>
    <tr>
      <td><strong>STL</strong></code></code></td>
      <td><code>PPC_NUM_THREADS=N</code></code></code></code></td>
      <td><code>export PPC_NUM_THREADS=4</code></code></code></code></td>
    </tr>
    <tr>
      <td><strong>ALL</strong></code></code></td>
      <td><code>OMP_NUM_THREADS=N</code>, <code>mpiexec -np P</code></code></code></code></td>
      <td><code>export OMP_NUM_THREADS=2</code><br><code>mpirun -np 2 ./program</code></code></code></code></td>
    </tr>
  </tbody>
</table>

### Входные тестовые данные

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center; width: 100%; font-family: Arial, sans-serif;">
  <thead>
    <tr style="background-color: #f2f2f2;">
      <th>Тип теста</th>
      <th>Размер матриц</th>
      <th>Плотность</th>
      <th>Количество ненулевых</th>
      <th>Тип заполнения</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>Диагональные матрицы</strong></code></code></td>
      <td>100,000 × 100,000</code></code></td>
      <td>0.0001%</code></code></td>
      <td>100,000</code></code></td>
      <td>Диагональное (1 элемент на строку)</code></code></td>
    </tr>
    <tr>
      <td><strong>Случайные матрицы</strong></code></code></td>
      <td>10,000 × 10,000</code></code></td>
      <td>0.1%</code></code></td>
      <td>~100 на строку (всего ~1,000,000)</code></code></td>
      <td>Случайное (равномерное распределение)</code></code></td>
    </tr>
  </tbody>
</table>

### Метрики производительности

`Speedup (ускорение)`: `Speedup = T_seq / T_parallel`

`T_seq` — время выполнения последовательной версии (task_run)

`T_parallel` — время выполнения параллельной версии

Efficiency (эффективность): `Efficiency = (Speedup / Workers) × 100%`

`Workers` — количество задействованных вычислительных единиц (потоков/процессов)

### Процедура измерений
- Каждый тест запускается многократно, фиксируется медианное время

- Прогрев выполняется перед замерами (инфраструктурой тестов)

- Окружение изолировано от фоновых процессов

## 4. Сводка корректности

### Метод верификации

Все параллельные реализации сравнивались с последовательной эталонной версией (SEQ) на одинаковых наборах входных данных. Проверка включала:

- Полное совпадение размеров результирующей матрицы

- Совпадение количества ненулевых элементов

- Совпадение значений с точностью `e = 1e-10`

- Валидность CRS-формата (`row_ptr` — неубывающая последовательность)

### Функциональные тесты (4 тестовых случая)

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>№</th> <th>Название</th> <th>Размеры</th> <th>Описание</th> <th>Ожидаемый nnz</th> </tr> </thead> <tbody> <tr> <td>1</code></code></td> <td>identity</code></code></td> <td>2×2 × 2×2</code></code></td> <td>A = [[1,0],[0,1]], B = [[2,0],[0,3]]</code></code></td> <td>2</code></code></td> </tr> <tr> <td>2</code></code></td> <td>simple_2x2</code></code></td> <td>2×2 × 2×2</code></code></td> <td>Разреженные матрицы с 3 ненулевыми в A</code></code></td> <td>3</code></code></td> </tr> <tr> <td>3</code></code></td> <td>zero_matrix</code></code></td> <td>2×2 × 2×2</code></code></td> <td>B — нулевая матрица</code></code></td> <td>0</code></code></td> </tr> <tr> <td>4</code></code></td> <td>sparse_3x3</code></code></td> <td>3×3 × 3×3</code></code></td> <td>Разреженные матрицы сложной структуры</code></code></td> <td>4</code></code></td> </tr> </tbody> </table>

### Результаты тестирования

- Все функциональные тесты пройдены для всех 5 реализаций

- Результаты параллельных версий совпадают с SEQ (погрешность < `1e-10`)

- Ограничения применимости: ни одна из реализаций не имеет принципиальных ограничений

## 5. Агрегированные результаты

### Тестирование на диагональных матрицах (100,000×100,000)

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Технология</th> <th>Конфигурация</th> <th>Time (task), s</th> <th>Speedup (vs SEQ)</th> <th>Efficiency</th> </tr> </thead> <tbody> <tr> <td><strong>SEQ</strong></code></code></td> <td>1 ранг</code></code></td> <td>0.0119372845</code></code></td> <td>1.00×</code></code> (baseline)</code></code></td> <td>100%</code></code></td> </tr> <tr> <td><strong>OMP</strong></code></code></td> <td>12 потоков</code></code></td> <td>0.0717499506</code></code></td> <td><strong style="color: #e74c3c;">0.16×</strong></code></code> (в 6.1× медленнее)</code></code></td> <td>1.4%</code></code></td> </tr> <tr> <td><strong>TBB</strong></code></code></td> <td>2 ранга</code></code></td> <td>0.0194890976</code></code></td> <td><strong style="color: #e67e22;">0.61×</strong></code></code> (в 1.6× медленнее)</code></code></td> <td>31%</code></code></td> </tr> <tr> <td><strong>STL</strong></code></code></td> <td>4 потока</code></code></td> <td>0.0326702604</code></code></td> <td><strong style="color: #e67e22;">0.37×</strong></code></code> (в 2.7× медленнее)</code></code></td> <td>9%</code></code></td> </tr> <tr> <td><strong style="color: #27ae60;">ALL (MPI+OMP)</strong></code></code></td> <td>2 процесса × 2 потока</code></code></td> <td><strong>0.0056650746</strong></code></code></td> <td><strong style="color: #27ae60;">2.14×</strong></code></code> </code></code></td> <td><strong style="color: #27ae60;">54%</strong></code></code></td> </tr> </tbody> </table>

### Тестирование на случайных матрицах (10,000×10,000, плотность 0.1%)

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Технология</th> <th>Конфигурация</th> <th>Time (task), s</th> <th>Speedup (vs SEQ)</th> <th>Efficiency</th> </tr> </thead> <tbody> <tr> <td><strong>SEQ</strong></code></code></td> <td>1 ранг</code></code></td> <td>0.260065</code></code></td> <td>1.00×</code></code> (baseline)</code></code></td> <td>100%</code></code></td> </tr> <tr> <td><strong>OMP</strong></code></code></td> <td>12 потоков</code></code></td> <td>0.077307</code></code></td> <td><strong style="color: #27ae60;">3.36×</strong></code></code> </code></code></td> <td>28%</code></code></td> </tr> <tr> <td><strong>TBB</strong></code></code></td> <td>12 потоков</code></code></td> <td>0.084945</code></code></td> <td><strong style="color: #27ae60;">3.06×</strong></code></code> </code></code></td> <td>26%</code></code></td> </tr> <tr> <td><strong>STL</strong></code></code></td> <td>4 потока</code></code></td> <td>0.084496</code></code></td> <td><strong style="color: #27ae60;">3.08×</strong></code></code> </code></code></td> <td><strong style="color: #27ae60;">77%</strong></code></code></td> </tr> </tbody> </table>

## 6. Интерпретация различий

### SEQ

Последовательная версия является эталоном корректности и базой для расчёта ускорения. Её детерминированность и простота верификации делают её идеальным baseline'ом для сравнения.

### OMP

#### Сильные стороны:

- Простота реализации — достаточно добавить несколько директив #pragma omp
- Динамическое планирование `(schedule(dynamic))` эффективно балансирует неравномерную нагрузку
- На случайных матрицах показал наилучшее ускорение (3.36×)

#### Слабые стороны:

- На диагональных матрицах проваливается из-за оверхэда (0.16×)
- Требует поддержки OpenMP компилятором

### TBB

#### Сильные стороны:

- Work-stealing планировщик адаптируется к неравномерной нагрузке

- На диагональных матрицах показал лучший результат среди однопоточных технологий (0.61×)

- Кроссплатформенность

#### Слабые стороны:

- Больший overhead на мелких задачах, чем у OMP

- На случайных матрицах уступает OMP (3.06× против 3.36×)

### STL
#### Сильные стороны:

- Наивысшая эффективность на случайных матрицах (77% при 4 потоках)

- Нет внешних зависимостей

- Полный контроль над потоками

#### Слабые стороны:

- Статическое разбиение неэффективно при неравномерной нагрузке

- Высокая сложность корректной реализации

- Плохое масштабирование на диагональных матрицах

### ALL (MPI + OpenMP)
### Сильные стороны:

- Единственная технология, достигшая ускорения >1× на диагональных матрицах (2.14×)

- Двухуровневый параллелизм компенсирует оверхэд даже на микроскопической нагрузке

- Возможность использования нескольких узлов

#### Слабые стороны:

- Высокая сложность реализации

- Дублирование матрицы B во всех процессах (расход памяти)

- Накладные расходы на MPI-коммуникации

### Сводная таблица сравнения

## 7. Репродуцируемость

### Команды сборки

```
# Сборка проекта
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release 
cmake --build build --parallel
```

### Команды запуска функциональных тестов

```
# SEQ
./build/bin/ppc_func_tests --gtest_filter="*tsyplakov*seq*"

# OMP (12 потоков)
export OMP_NUM_THREADS=12
./build/bin/ppc_func_tests --gtest_filter="*tsyplakov*omp*"

# TBB (12 потоков)
export PPC_NUM_THREADS=12
./build/bin/ppc_func_tests --gtest_filter="*tsyplakov*tbb*"

# STL (4 потока)
export PPC_NUM_THREADS=4
./build/bin/ppc_func_tests --gtest_filter="*tsyplakov*stl*"

# ALL (2 процесса × 2 потока)
export OMP_NUM_THREADS=2
mpirun -np 2 ./build/bin/ppc_func_tests --gtest_filter="*tsyplakov*all*"
```

### Команда запуска производительных тестов

```
# SEQ (baseline)
./build/bin/ppc_perf_tests --gtest_filter="*tsyplakov*seq*"

# OMP (12 потоков)
export OMP_NUM_THREADS=12
./build/bin/ppc_perf_tests --gtest_filter="*tsyplakov*omp*"

# TBB (12 потоков)
export PPC_NUM_THREADS=12
./build/bin/ppc_perf_tests --gtest_filter="*tsyplakov*tbb*"

# STL (4 потока)
export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*tsyplakov*stl*"

# ALL (2 процесса × 2 потока)
export OMP_NUM_THREADS=2
mpirun -np 2 ./build/bin/ppc_perf_tests --gtest_filter="*tsyplakov*all*"
```

## 8. Заключение

В рамках данной работы успешно реализованы и протестированы пять вариантов умножения разреженных матриц. Все реализации корректны и проходят функциональные тесты.

### Основные выводы

1. ```OMP``` - Лучшее соотношение производительности к сложности реализации.
2. `TBB` - Кроссплатформенные задачи, где важна переносимость.
3. `STL` - Максимальный контроль и минимальный `overhead` (до 4 потоков).
4. `ALL` - Для суперкомпьютеров, кластеров, очень крупных задач.

### Ключевые результаты

1. На случайных матрицах (реалистичная нагрузка) все три потоковые технологии показали ускорение >3×:

- OpenMP: 3.36× (12 потоков)

- TBB: 3.06× (12 потоков)

- STL: 3.08× (4 потока, эффективность 77%)

2. На диагональных матрицах:

- OpenMP провалился (0.16×)

- TBB и STL показали замедление (0.61× и 0.37×)

- ALL выиграл (2.14×) — единственный, кто справился с микроскопической нагрузкой

3. Эффективность STL на 4 потоках (77%) оказалась лучшей среди всех технологий, что подтверждает ценность ручного управления при правильном подходе.

