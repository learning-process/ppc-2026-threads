# Умножение разреженных матриц в формате CRS (Compressed Row Storage) — ALL (MPI + OpenMP)
- Студент: Цыплаков Кирилл
- Технология: ALL (MPI + OpenMP)
- Вариант: 4

## 1. Контекст
Гибридная технология ALL сочетает два уровня параллелизма:

- Межпроцессный уровень (MPI) — распределение работы между независимыми процессами. Каждый процесс имеет свою собственную память и может работать на отдельном ядре/узле.

- Внутрипроцессный уровень (OpenMP) — распараллеливание внутри каждого процесса с использованием общей памяти.

Такая двухуровневая модель особенно эффективна для:

- Крупных задач, не помещающихся в память одного узла

- Систем с гибридной архитектурой (многоядерные узлы + кластеризация)

- Задач с естественным разделением данных, где каждый процесс может обрабатывать свой блок независимо

В данной реализации ALL применяется для умножения разреженных матриц: MPI распределяет строки матрицы A между процессами, а внутри каждого процесса OpenMP распараллеливает обработку строк.

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

### Особенности ALL-версии
- Двухуровневый параллелизм: MPI между процессами + OpenMP внутри процесса

- Каждый процесс получает свою порцию строк матрицы A

- Матрица B реплицируется (копируется) во все процессы (только чтение)

- Внутри каждого процесса используется OpenMP для распараллеливания строк

- Сбор результатов от всех процессов в корневой процесс (rank 0)

## 3. Базовый алгоритм

### Принцип работы
Гибридный алгоритм комбинирует:

- MPI уровень: Распределение строк матрицы A между процессами

- OpenMP уровень: Внутри каждого процесса — параллельная обработка полученных строк

### Пошаговое описание
#### Шаг 1: Инициализация MPI
Определяются ранг процесса `(rank)` и общее количество процессов `(size)`.

#### Шаг 2: Распределение строк
Матрица `A` разбивается на `size` непрерывных блоков строк. Каждый процесс получает свой блок.

#### Шаг 3: Репликация матрицы B
Матрица `B` отправляется во все процессы (или уже присутствует, если входные данные реплицированы).

#### Шаг 4: Локальное умножение (OpenMP)
Каждый процесс выполняет умножение своей порции строк матрицы `A` на матрицу `B` с использованием OpenMP.

#### Шаг 5: Сбор результатов
Все процессы отправляют свои результаты (части матрицы C) в процесс с рангом 0.

#### Шаг 6: Формирование выходной матрицы
Процесс 0 объединяет полученные части в итоговую матрицу C.

### Асимптотическая сложность

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center; width: 100%;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Метрика</th> <th>Формула</th> <th>Примечание</th> </tr> </thead> <tbody> <tr> <td><strong>Время (параллельное)</strong></code></code></td> <td><code>T_parallel ≈ T_compute / (P_mpi × P_omp) + T_comm + T_overhead</code></code></td> <td>T_comm — время на MPI коммуникации (распределение + сбор)</code></code></td> </tr> <tr> <td><strong>Время (коммуникации)</strong></code></code></td> <td><code>T_comm ∝ P_mpi × (размер сообщения)</code></code></td> <td>Линейно зависит от числа процессов</code></code></td> </tr> <tr> <td><strong>Время (худший случай)</strong></code></code></td> <td><code>O(n³ / (P_mpi × P_omp)) + O(P_mpi × n²)</code></code></td> <td>Первое слагаемое — вычисления, второе — коммуникации</code></code></td> </tr> </tbody> </table>

### Критерий корректности

Результат умножения должен совпадать с результатом умножения плотных матриц, преобразованных из CRS-формата. Допустимая погрешность: `e = 1e-12`.

## 4. Схема распараллеливания

### MPI уровень (межпроцессный)

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Параметр</th> <th>Значение</th> <th>Обоснование</th> </tr> </thead> <tbody> <tr> <td><strong>Распределение строк</strong></code></code></td> <td>Статическое (непрерывные блоки)</code></code></td> <td>Каждый процесс получает последовательный диапазон строк матрицы A</code></code></td> </tr> <tr> <td><strong>Коммуникация B</strong></code></code></td> <td><code>MPI_Bcast</code> (широковещательная рассылка)</code></code></td> <td>Матрица B реплицируется во все процессы для локального доступа</code></code></td> </tr> <tr> <td><strong>Сбор результатов</strong></code></code></td> <td><code>MPI_Gather</code> / <code>MPI_Send</code> + <code>MPI_Recv</code></code></code></td> <td>Каждый процесс отправляет свою часть C в корневой процесс (rank 0)</code></code></td> </tr> <tr> <td><strong>Синхронизация</strong></code></code></td> <td><code>MPI_Barrier</code> (при необходимости)</code></code></td> <td>Гарантирует готовность данных перед началом вычислений</code></code></td> </tr> </tbody> </table>

### OpenMP уровень (внутрипроцессный)

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Параметр</th> <th>Значение</th> <th>Обоснование</th> </tr> </thead> <tbody> <tr> <td><strong>Распараллеливаемая область</strong></code></code></td> <td>Цикл по локальным строкам матрицы A</code></code></td> <td>Вычисления для каждой строки полностью независимы</code></code></td> </tr> <tr> <td><strong>Директива</strong></code></code></td> <td><code>#pragma omp parallel for schedule(dynamic)</code></code></code></td> <td>Параллельное выполнение итераций с динамической балансировкой</code></code></td> </tr> <tr> <td><strong>Распределение итераций</strong></code></code></td> <td><code>schedule(dynamic)</code></code></code></td> <td>Компенсирует возможную неравномерность плотности строк</code></code></td> </tr> <tr> <td><strong>Атрибуты переменных</strong></code></code></td> <td><code>default(none) shared(a, b, row_values, row_cols)</code></code></code></td> <td>Явное указание атрибутов повышает безопасность</code></code></td> </tr> <tr> <td><strong>Синхронизация</strong></code></code></td> <td>Неявный барьер в конце <code>parallel for</code></code></code></td> <td>Ожидание завершения всех потоков перед продолжением</code></code></td> </tr> </tbody> </table>

### Переменные и их доступ

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: left; width: 100%;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Переменная</th> <th>Уровень</th> <th>Доступ</th> <th>Обоснование</th> </tr> </thead> <tbody> <tr> <td><code>A</code></code></td> <td>MPI (распределена)</code></code></td> <td>Каждый процесс — свой блок</code></code></td> <td>Строки разделены между процессами</code></code></td> </tr> <tr> <td><code>B</code></code></td> <td>MPI (реплицирована)</code></code></td> <td>Только чтение во всех процессах</code></code></td> <td>Репликация для локального доступа</code></code></td> </tr> <tr> <td><code>local_C</code></code></td> <td>Локальная для процесса</code></code></td> <td>Каждый процесс пишет в свою часть</code></code></td> <td>Гонок между процессами нет</code></code></td> </tr> <tr> <td><code>row_values</code></code> (OpenMP)</code></code></td> <td>Внутри процесса</code></code></td> <td>Разные потоки — разные строки</code></code></td> <td>Гонок внутри процесса нет</code></code></td> </tr> </tbody> </table>

## 5. Детали реализации

### Файлы
`all/include/ops_all.hpp` — заголовочный файл

`all/src/ops_all.cpp` — реализация

### Ключевые фрагменты кода

Основная логика RunImpl (псевдокод на основе типовой реализации):

```
bool TsyplakovKTestTaskALL::RunImpl() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    const auto &a = GetInput().a;
    const auto &b = GetInput().b;
    
    // Репликация матрицы B во все процессы
    SparseMatrixCRS local_b = b;
    MPI_Bcast(&local_b, ...);
    
    // Распределение строк A между процессами
    int rows_per_proc = a.rows / size;
    int start_row = rank * rows_per_proc;
    int end_row = (rank == size - 1) ? a.rows : start_row + rows_per_proc;
    
    // Локальное умножение с OpenMP
    SparseMatrixCRS local_c(a.rows, b.cols);
    
    #pragma omp parallel for schedule(dynamic)
    for (int i = start_row; i < end_row; ++i) {
        ComputeRow(a, b, i, local_c);
    }
    
    // Сбор результатов в процесс 0
    if (rank == 0) {
        // Объединение всех local_c в итоговую C
        GatherResults(local_c, ...);
    } else {
        MPI_Send(local_c, ..., 0, ...);
    }
    
    return true;
}
```
Функция ComputeRow (аналогична STL/TBB версиям):

```
void ComputeRow(const SparseMatrixCRS &a, const SparseMatrixCRS &b,
                int row, SparseMatrixCRS &c) {
    std::unordered_map<int, double> acc;
    
    for (int idx_a = a.row_ptr[row]; idx_a < a.row_ptr[row + 1]; ++idx_a) {
        const int k = a.col_index[idx_a];
        const double val_a = a.values[idx_a];
        
        for (int idx_b = b.row_ptr[k]; idx_b < b.row_ptr[k + 1]; ++idx_b) {
            const int j = b.col_index[idx_b];
            acc[j] += val_a * b.values[idx_b];
        }
    }
    
    for (const auto &[col, val] : acc) {
        if (std::fabs(val) > 1e-12) {
            c.values.push_back(val);
            c.col_index.push_back(col);
        }
    }
}
```

### Особенности реализации
1. Двухуровневый параллелизм:

- MPI распределяет строки между процессами

- OpenMP распределяет строки внутри процесса между потоками

- Общая степень параллелизма: P_mpi × P_omp

2. Репликация матрицы B:

- B копируется во все процессы для локального доступа

- Избегает дополнительных коммуникаций при вычислениях

3. Динамическое распределение в OpenMP:

- schedule(dynamic) для балансировки нагрузки

- Важно при нерегулярной плотности строк

4. Сбор результатов:

- Каждый процесс отправляет свою часть матрицы C в процесс 0

- Процесс 0 объединяет результаты

## 6. Проверка корректности

### Метод верификации

Корректность ALL-версии проверяется сравнением с SEQ:

- Все функциональные тесты прогоняются на обеих реализациях

- Результаты сравниваются поэлементно с погрешностью `ε = 1e-10`

- Дополнительно проверяется согласованность результатов между разными MPI процессами

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

Тесты производительности запускались на матрицах размера 100,000×100,000 с диагональным заполнением. Плотность ненулевых элементов составляет 0.0001% (1 ненулевой элемент на строку).

### Полная сводная таблица всех технологий

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; text-align: center;"> <thead> <tr style="background-color: #f2f2f2;"> <th>Технология</th> <th>Mode</th> <th>Time, s</th> <th>Speedup (vs SEQ task)</th> </tr> </thead> <tbody> <tr> <td><strong>SEQ</strong></code></code></td> <td>task_run</code></code></td> <td>0.0121132356</code></code></td> <td><strong>1.00×</strong></code></code> (baseline)</code></code></td> </tr> <tr> <td><strong>SEQ</strong></code></code></td> <td>pipeline</code></code></td> <td>0.0135521440</code></code></td> <td>0.89×</code></code></td> </tr> <tr> <td><strong>OMP</strong></code></code></td> <td>task_run</code></code></td> <td>0.0691685780</code></code></td> <td><strong style="color: #e74c3c;">0.18×</strong> (в 5.7× медленнее)</code></code></td> </tr> <tr> <td><strong>OMP</strong></code></code></td> <td>pipeline</code></code></td> <td>0.0437379110</code></code></td> <td>0.28×</code></code></td> </tr> <tr> <td><strong>TBB</strong></code></code></td> <td>task_run</code></code></td> <td>0.0188774118</code></code></td> <td>0.64×</code></code></td> </tr> <tr> <td><strong>TBB</strong></code></code></td> <td>pipeline</code></code></td> <td>0.0200403730</code></code></td> <td>0.60×</code></code></td> </tr> <tr> <td><strong>STL</strong></code></code></td> <td>task_run</code></code></td> <td>0.0326702604</code></code></td> <td>0.37×</code></code></td> </tr> <tr> <td><strong>STL</strong></code></code></td> <td>pipeline</code></code></td> <td>0.0207790404</code></code></td> <td>0.58×</code></code></td> </tr> <tr> <td><strong style="color: #27ae60;">ALL (MPI+OMP)</strong></code></code></td> <td>pipeline</code></code></td> <td>0.0064839572</code></code></td> <td><strong style="color: #27ae60;">1.87×</strong></code></code> </code></code></td> </tr> <tr> <td><strong style="color: #27ae60;">ALL (MPI+OMP)</strong></code></code></td> <td>task_run</code></code></td> <td><strong>0.0056650746</strong></code></code></td> <td><strong style="color: #27ae60;">2.14×</strong></code></code> </code></code></td> </tr> </tbody> </table>

### Анализ производительности

Ключевое наблюдение: ALL-версия — единственная технология, показавшая реальное ускорение на диагональных матрицах.

1. Двухуровневый параллелизм компенсирует оверхэд:

- MPI распределяет строки между процессами

- OpenMP распределяет строки внутри процесса

- Даже при микроскопической нагрузке (1 операция на строку) общий оверхэд распределяется между worker'ами

2. Эффект разделения:

- 2 MPI процесса × 2 OpenMP потока = 4 worker'а

- Каждый worker обрабатывает 25,000 строк

- Оверхэд на worker становится меньше

3. Репликация `B` не создаёт проблем:

- `B` маленькая (1 ненулевой на строку)

- Коммуникационные затраты минимальны

4. Сбор результатов эффективен:

- Каждый процесс отправляет только свою часть `C`

- При малом количестве процессов оверхэд незначителен

## 9. Выводы

ALL — единственная технология, достигшая ускорения >1× на диагональных матрицах.

```Speedup: 2.14× (task_run) и 1.87× (pipeline)```

Это лучший результат среди всех пяти реализаций.

Почему ALL победила:

1. Двухуровневый параллелизм (MPI + OpenMP) лучше распределяет даже микроскопическую нагрузку

2. Оверхэд распределяется между несколькими процессами и потоками

3. Коммуникационные затраты минимальны для диагональных матриц

### Границы применимости
**ALL-версия будет эффективна, когда**:

- Задача достаточно велика (чтобы окупить MPI коммуникации)

- Доступно несколько MPI процессов (кластер или многоядерная система)

- Требуется максимальная производительность на диагональных/разреженных матрицах

- Память одного узла недостаточна для всей задачи

**ALL-версия менее эффективна, когда**:

- Задача мала (оверхэд MPI доминирует)

- Матрицы очень плотные (коммуникации могут стать узким местом)

- Нет поддержки MPI в окружении

### Итог
ALL-реализация показала, что гибридный подход (MPI + OpenMP) является наиболее мощным инструментом для параллельного умножения разреженных матриц, особенно на диагональных/разреженных данных. Даже при микроскопической вычислительной нагрузке ALL смогла достичь ускорения 2.14×, в то время как все остальные технологии показали замедление.

Это подтверждает, что для максимальной производительности на современных гибридных системах (многоядерные узлы + кластеризация) необходимо использовать двухуровневый параллелизм.