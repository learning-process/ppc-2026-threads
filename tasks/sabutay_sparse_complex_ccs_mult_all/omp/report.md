# Умножение разреженных комплексных матриц (CCS) — OMP

- Student: Иманов Сабутай Ширзад оглы, группа 3823Б1ПР5
- Technology: OMP
- Variant: 7

## 1. Контекст

В этой версии умножение матриц в формате CCS выполняется с помощью **OpenMP**: несколько потоков обрабатывают разные
столбцы результата одновременно. Для сравнения времени используется эталон SEQ: $T_{seq} = 0.00001390$ с при одном
потоке (`task_run`).

## 2. Постановка задачи

Нужно вычислить произведение $C = A \cdot B$ для двух разреженных **комплексных** матриц в формате **CCS**.

- **Вход:** `InType = std::tuple<CCS, CCS>`.
- **Выход:** `OutType = CCS`.
- **Размеры:** $\text{cols}(A) = \text{rows}(B)$ — проверяется в `ValidationImpl` класса
`SabutaySparseComplexCcsMultOmpFix`.
- **Структура CCS:** корректность массивов `col_start`, `row_index`, `nz` и индексов строк.
- **Слияние:** отбрасывание малых значений ($|z| \le 10^{-14}$).
- **Эталон времени:** $T_{seq} = 0.00001390$ с (SEQ, 1 поток, `task_run`).

## 3. Базовый алгоритм

Столбцовое умножение разреженных матриц:

1. для столбца $j$ матрицы $B$ — собрать вклады (`BuildColumnFromRight`);
2. отсортировать буфер по номеру строки (`SortBufferByRow`);
3. слить повторяющиеся строки (`CoalesceBufferToColumn` / сборка в `out`).

В варианте **OMP** шаги 1–3 для разных $j$ выполняются **параллельно** в `BuildProductMatrixOmp`
(`all/src/ops_all.cpp`). После параллельной фазы столбцы **последовательно** копируются в итоговую CCS.

## 4. Схема распараллеливания

- Создаётся команда потоков OpenMP (`#pragma omp parallel`).
- Столбцы результата делятся между потоками (`#pragma omp for schedule(static)`).
- У каждого потока свой временный буфер — потоки не мешают друг другу при расчёте столбца.
- После параллельной части один поток (главный) последовательно собирает все столбцы в итоговую матрицу `out`.

Число потоков задаётся `PPC_NUM_THREADS` (курс дублирует в `OMP_NUM_THREADS`).

**Области данных OpenMP:**

| Переменная / объект | Атрибут | Пояснение |
| ------------------- | ------- | --------- |
| `left`, `right`, `columns` | `shared` | общие входы и массив столбцов |
| `buffer` | фактически **private** | объявлен внутри тела `parallel` |
| `j` | private | индекс цикла `for` |

`reduction` **не нужен**: каждый столбец $j$ пишет только в `columns[j]`. В конце области `parallel` есть **неявный
барьер**; затем идёт последовательная склейка в `out`.

`schedule(static)` — столбцы делятся заранее поровну; подходит, когда работа по столбцам примерно одинакова.

## 5. Детали реализации

Файлы: `omp/src/ops_omp.cpp`, общая функция умножения — `BuildProductMatrixOmp` в `all/src/ops_all.cpp`.

Класс OMP вызывает умножение с типом `kOMP`:

```55:63:tasks/sabutay_sparse_complex_ccs_mult_all/omp/src/ops_omp.cpp
void SabutaySparseComplexCcsMultOmpFix::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out) {
  SabutaySparseComplexCcsMultAll::BuildProductMatrix(left, right, out, ppc::task::TypeOfTask::kOMP);
}

bool SabutaySparseComplexCcsMultOmpFix::RunImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  BuildProductMatrix(left, right, GetOutput());
  return true;
}
```

### Листинг: параллельный цикл по столбцам

```163:177:tasks/sabutay_sparse_complex_ccs_mult_all/all/src/ops_all.cpp
#pragma omp parallel default(none) shared(left, right, columns)
  {
    std::vector<std::pair<int, std::complex<double>>> buffer;
    buffer.reserve(128U);

#pragma omp for schedule(static)
    for (int j = 0; j < right.col_count; ++j) {
      auto &column = columns[static_cast<std::size_t>(j)];
      BuildColumnFromRight(left, right, j, buffer);
      if (!buffer.empty()) {
        SortBufferByRow(buffer);
      }
      CoalesceBufferToColumn(buffer, column, kDropMagnitude);
    }
  }
```

_Пояснение:_ каждый поток считает свои столбцы и кладёт результат во временный массив `columns`. Затем в обычном цикле
столбцы копируются в итоговую CCS.

## 6. Проверка корректности

Функциональные тесты из `tests/functional/main.cpp`: результат сравнивается с плотным умножением, допуск $10^{-12}$.
Класс — `SabutaySparseComplexCcsMultOmpFix`.

**Результаты func-тестов:** 3/3 **PASSED**.

| `case_id` | Суффикс GTest | Результат |
| :-------: | ------------- | :-------: |
| 0 | `..._omp_enabled_0` | PASSED |
| 1 | `..._omp_enabled_1` | PASSED |
| 2 | `..._omp_enabled_2` | PASSED |

```powershell
$env:PPC_NUM_PROC = "1"
mpiexec -n 1 build\bin\ppc_func_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*omp_enabled*"
```

## 7. Экспериментальная среда

- **CPU:** Intel Core Ultra 5 125H, 14 ядер / 18 потоков.
- **RAM:** 32 ГБ, **ОС:** Windows 11, **сборка:** Release (MSVC).
- **Размер задачи (perf):** случайные матрицы $80 \times 90$ и $90 \times 70$ (`tests/performance/main.cpp`).
- **Прогонов внутри замера:** 5.
- **Команда:**

```powershell
$env:PPC_NUM_PROC = "1"
$env:PPC_NUM_THREADS = "4"   # для каждого из 1, 2, 4, 8
build\bin\ppc_perf_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*omp_enabled*"
```

## 8. Результаты

Каркас perf запускает **`task_run`** и **`pipeline`**. Для $S$ и $E$ ниже: $T_{seq} = 0.00001390$ с, $S = T_{seq} / T$,
$E = S / p \cdot 100\%$.

### `task_run`

| $p$ | Время (с) | $S$ | $E$ |
| :-: | :-------: | :-: | :-: |
| 1 | 0.00002268 | 0.61 | 61% |
| 2 | 0.00008498 | 0.16 | 8% |
| 4 | 0.00007856 | 0.18 | 5% |
| 8 | 0.00011186 | 0.12 | 2% |

_Комментарий:_ на тестовых матрицах из perf-набора столбцов немного (около 70), а запуск OpenMP стоит дорого. Поэтому
параллельная версия **не быстрее** SEQ — это нормально для маленькой задачи. На больших матрицах выигрыш обычно
заметнее.

### `pipeline`

| $p$ | Время (с) |
| :-: | :-------: |
| 1 | 0.00044970 |
| 2 | 0.00076444 |
| 4 | 0.00177710 |
| 8 | 0.00148146 |

## 9. Выводы

OMP распределяет столбцы между потоками — это понятная схема для CCS. В perf-тесте курса ускорения нет из-за малого
размера данных. Реализация совпадает с эталоном SEQ по результату.
