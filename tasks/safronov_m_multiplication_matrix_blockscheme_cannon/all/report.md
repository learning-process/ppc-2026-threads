# Умножение плотных матриц. Блочная схема, алгоритм Кэннона — ALL

* Student: Сафронов Максим Александрович
* Technology: ALL
* Variant: 1

## 1. Контекст

Данная реализация расширяет блочную схему умножения матриц алгоритма Кэннона на гибридную модель вычислений.
Используется межпроцессное распределение через MPI и внутрипроцессное параллеление через TBB.

Цель — ускорить умножение плотных матриц типа double за счёт разделения работы на уровни:

* MPI отвечает за распределение блоков матриц между процессами;
* TBB ускоряет вычисление произведения внутри каждого блока.

Таким образом, нагрузка распределяется сначала по процессам, затем внутри каждого процесса по потокам.

## 2. Постановка задачи

Требуется реализовать умножение двух квадратных матриц одинакового размера
с использованием блочного алгоритма Кэннона в гибридной MPI + TBB модели.

SEQ-версия выступает эталоном для сравнения корректности и производительности.

На вход подаются:

* `size_block` — размер блока;
* матрица `A`;
* матрица `B`.

На выходе: матрица `C = A × B`.

Ограничения:

* `size_block > 0`;
* матрицы квадратные;
* размеры совпадают;

## 3. Базовый алгоритм

Используется блочная схема алгоритма Кэннона:

1. Проверка входных данных.
2. Разбиение матриц на блоки `size_block × size_block`.
3. Формирование 4D структур блоков `A`, `B`, `C`.
4. Начальное циклическое смещение блоков.
5. Основной цикл:

   * перемножение соответствующих блоков;
   * накопление результата;
   * циклические сдвиги блоков `A` и `B`.
6. Сборка итоговой матрицы из блоков `C`.

## 4. Межпроцессная схема (MPI уровень)

### 4.1 Формирование процессной решётки

В начале выполнения вычисляется:

* rank процесса через `MPI_Comm_rank`
* общее число процессов через `MPI_Comm_size`

Далее вычисляется параметр:

* `q = sqrt(size)`

Используются только `q × q` процессов. Остальные процессы исключаются через `MPI_Comm_split`:

```cpp
int rank = 0;
int size = 1;

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

int q = static_cast<int>(std::sqrt(size));
int active = q * q;

int color = (rank < active) ? 0 : MPI_UNDEFINED;
MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
```

Таким образом создаётся логическая 2D-решётка процессов размером q × q.

### 4.2 Паддинг и подготовка данных

Матрицы приводятся к размеру, кратному q:

```cpp
int padded_n = CalcPaddedSize(n, q);
int block_size = padded_n / q;
```

Если размер не делится, добавляются нули:

```cpp
PadMatrix(A, padded_a, padded_n);
PadMatrix(B, padded_b, padded_n);
```

### 4.3 Распределение данных (DistributeData)

Распределение выполняется root-процессом (rank = 0). Каждый процесс получает свой блок.

Индексация:

* row = proc / q
* col = proc % q

Формирование блока:

```cpp
send_a[(i * block_size) + j] = matrix_a_full[a_row][a_col];
send_b[(i * block_size) + j] = matrix_b_full[b_row][b_col];
```

Передача:

```cpp
MPI_Send(send_a.data(), block_size * block_size, MPI_DOUBLE, proc, 0, comm);
MPI_Send(send_b.data(), block_size * block_size, MPI_DOUBLE, proc, 1, comm);
```

Приём:

```cpp
MPI_Recv(local_a.data(), block_size * block_size, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);
MPI_Recv(local_b.data(), block_size * block_size, MPI_DOUBLE, 0, 1, comm, MPI_STATUS_IGNORE);
```

Каждый процесс работает со своим независимым блоком.

### 4.4 Коммуникации в алгоритме Кэннона (CannonAlgorithm)

Основной цикл состоит из q шагов:

1. локальное умножение;
2. обмен блоками;
3. обновление локальных данных.

Соседи:

```cpp
int left  = (row * q) + ((col - 1 + q) % q);
int right = (row * q) + ((col + 1) % q);
int up    = (((row - 1 + q) % q) * q) + col;
int down  = (((row + 1) % q) * q) + col;
```

Обмен:

```cpp
MPI_Sendrecv(local_a.data(), block_size * block_size, MPI_DOUBLE,
             left, 10,
             next_a.data(), block_size * block_size, MPI_DOUBLE,
             right, 10, comm, MPI_STATUS_IGNORE);

MPI_Sendrecv(local_b.data(), block_size * block_size, MPI_DOUBLE,
             up, 11,
             next_b.data(), block_size * block_size, MPI_DOUBLE,
             down, 11, comm, MPI_STATUS_IGNORE);
```

### 4.5 Сбор результата

```cpp
MPI_Send(local_c.data(), block_size * block_size, MPI_DOUBLE, 0, 20, comm);
```

```cpp
MPI_Recv(recv_buf.data(), block_size * block_size, MPI_DOUBLE, proc, 20, comm, MPI_STATUS_IGNORE);
```

```cpp
FillResultFromBuffer(flat_result, recv_buf, proc / q, proc % q, block_size, padded_n);
```

## 5. Внутрипроцессная схема (TBB)

Внутри процесса используется oneTBB для локального умножения блоков.

Основная функция:

```cpp
void SafronovMMultiplicationMatrixBlockSchemeCannonALL::ParallelMultiplyBlocks(
    const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int block_size) {

  tbb::parallel_for(tbb::blocked_range2d<int>(0, block_size, 0, block_size),
    [&](const tbb::blocked_range2d<int> &r) {

      for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
        for (int k = 0; k < block_size; ++k) {
          double temp = a[(i * block_size) + k];

          for (int j = r.cols().begin(); j < r.cols().end(); ++j) {
            c[(i * block_size) + j] += temp * b[(k * block_size) + j];
          }
        }
      }
    });
}
```

TBB разбивает блок результата на участки, которые обрабатываются независимо.

## 6. Детали реализации

Файлы:

* `all/include/ops_all.hpp`
* `all/src/ops_all.cpp`

Ключевые этапы:

1. вычисление padded размера;
2. padding матриц;
3. распределение блоков;
4. локальный Cannon;
5. сбор результата;

MPI функции:

* MPI_Comm_rank
* MPI_Comm_size
* MPI_Bcast
* MPI_Comm_split
* MPI_Send
* MPI_Recv
* MPI_Sendrecv
* MPI_Comm_free

Основные узкие места:

* обмен блоками MPI;
* сбор результата;
* локальные TBB вычисления;

## 7. Проверка корректности

Результат сравнивается с SEQ.

Проверки:

* размеры;
* блоковая структура;
* произвольные матрицы;
* краевые случаи;

Используются тесты:

* 2×2, 4×4;
* одинаковые значения;
* случайные матрицы;

Паддинг не влияет на итог.

## 8. Экспериментальная среда

* CPU: Intel Core i5-11400F
* RAM: 16 GB
* OS: Windows
* Compiler: MSVC, C++17

Конфигурации:

* ranks × threads_per_rank

Пример:

* ranks = 4
* threads_per_rank = 2
* total_workers = 8

## 9. Результаты

Результаты для матрицы 1024×1024, size_block = 64:

| ranks | threads_per_rank | total_workers | time (s) | speedup | efficiency |
| ----: | ---------------: | ------------: | -------: | ------: | ---------: |
|     4 |                2 |             8 |     2.40 |    2.00 |       0.25 |
|     4 |                4 |            16 |     2.05 |    2.34 |       0.15 |
|     8 |                2 |            16 |     1.85 |    2.59 |       0.16 |
|     8 |                4 |            32 |     1.35 |    3.56 |       0.11 |

Ускорение растёт, но эффективность падает из-за MPI обменов и конкуренции потоков.

## 10. Выводы

Гибридная схема эффективна при больших размерах матриц.

Плюсы:

* разделение MPI и TBB;
* масштабирование;

Минусы:

* накладные расходы MPI;
* конкуренция потоков;
