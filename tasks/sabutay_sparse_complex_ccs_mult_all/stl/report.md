# Умножение разреженных комплексных матриц (CCS) — STL

- Student: Иманов Сабутай Ширзад оглы, группа 3823Б1ПР5
- Technology: STL
- Variant: 7

## 1. Контекст

Вариант **STL** в этой задаче — это класс `SabutaySparseComplexCcsMultSTL` в каркасе курса. Само умножение матриц
выполняется **последовательно**, тем же кодом, что и SEQ. Отдельные потоки `std::thread` показаны в технологии ALL
(после умножения), а не в STL-файле.

Эталон: $T_{seq} = 0.00001390$ с.

## 2. Постановка задачи

Нужно вычислить произведение $C = A \cdot B$ для двух разреженных **комплексных** матриц в формате **CCS**.

- **Вход:** `InType = std::tuple<CCS, CCS>`.
- **Выход:** `OutType = CCS`.
- **Проверка:** `ValidationImpl` класса `SabutaySparseComplexCcsMultSTL`.
- **Условие размеров:** $\text{cols}(A) = \text{rows}(B)$.
- **Порог слияния:** $|z| \le 10^{-14}$ — элемент не записывается в результат.
- **Эталон времени:** $T_{seq} = 0.00001390$ с.

## 3. Базовый алгоритм

В STL вызывается последовательная функция `SpmmAbc` (ветка `kSTL` в `BuildProductMatrix`):

1. цикл по столбцам $j$ правой матрицы $B$;
2. для каждого $j$ — `BuildColumnFromRight`, сортировка, слияние в столбец $C$;
3. обновление массива `col_start` для результата.

Параллельных диапазонов и `std::thread` в `stl/src/ops_stl.cpp` **нет** — умножение полностью однопоточное.

## 4. Схема распараллеливания

В **умножении** STL параллелизма нет: вызывается `SpmmAbc`, как у SEQ.

Потоки `std::thread` в этой задаче используются в классе **ALL** (`RunImpl`: создаётся `num_threads` потоков, каждый
увеличивает счётчик, затем `join`). В `stl/src/ops_stl.cpp` ручного разбиения диапазонов и `join` для SpMM **нет** — это
важно указать в отчёте честно.

`PPC_NUM_THREADS` на время умножения **не влияет**; в perf-тесте меняется только для единой методики курса.

## 5. Детали реализации

```45:63:tasks/sabutay_sparse_complex_ccs_mult_all/stl/src/ops_stl.cpp
void SabutaySparseComplexCcsMultSTL::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out) {
  SabutaySparseComplexCcsMultAll::BuildProductMatrix(left, right, out, ppc::task::TypeOfTask::kSTL);
}

bool SabutaySparseComplexCcsMultSTL::RunImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  BuildProductMatrix(left, right, GetOutput());
  return true;
}
```

В общем файле SEQ и STL идут в одну ветку:

```257:261:tasks/sabutay_sparse_complex_ccs_mult_all/all/src/ops_all.cpp
    case ppc::task::TypeOfTask::kSEQ:
    case ppc::task::TypeOfTask::kSTL:
      SpmmAbc(left, right, out, kDropMagnitude);
      return;
```

## 6. Проверка корректности

`tests/functional/main.cpp`, класс `SabutaySparseComplexCcsMultSTL`.

**Результаты func-тестов:** 3/3 **PASSED**.

| `case_id` | Суффикс GTest | Результат |
| :-------: | ------------- | :-------: |
| 0 | `..._stl_enabled_0` | PASSED |
| 1 | `..._stl_enabled_1` | PASSED |
| 2 | `..._stl_enabled_2` | PASSED |

```powershell
$env:PPC_NUM_PROC = "1"
mpiexec -n 1 build\bin\ppc_func_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*stl_enabled*"
```

## 7. Экспериментальная среда

- **CPU:** Intel Core Ultra 5 125H, 14 ядер / 18 потоков.
- **RAM:** 32 ГБ, **ОС:** Windows 11, **сборка:** Release (MSVC).
- **Размер задачи (perf):** матрицы $80 \times 90$ и $90 \times 70$.
- **Прогонов внутри замера:** 5.

Команда:

```powershell
$env:PPC_NUM_PROC = "1"
$env:PPC_NUM_THREADS = "4"
build\bin\ppc_perf_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*stl_enabled*"
```

## 8. Результаты

Каркас perf запускает **`task_run`** и **`pipeline`**. $T_{seq} = 0.00001390$ с.

### `task_run`

Значения $S > 1$ — погрешность замера, а не реальное ускорение (код однопоточный).

| $p$ | Время (с) | $S$ | $E$ |
| :-: | :-------: | :-: | :-: |
| 1 | 0.00001266 | 1.10 | 110% |
| 2 | 0.00001302 | 1.07 | 53% |
| 4 | 0.00001258 | 1.11 | 28% |
| 8 | 0.00001262 | 1.10 | 14% |

_Комментарий:_ время около $10^{-5}$ с и почти не зависит от $p$ — так и должно быть, потому что код однопоточный.
Небольшие отличия от SEQ — погрешность замера, а не ускорение.

### `pipeline`

| $p$ | Время (с) |
| :-: | :-------: |
| 1 | 0.00003026 |
| 2 | 0.00010842 |
| 4 | 0.00002598 |
| 8 | 0.00002726 |

## 9. Выводы

STL в проекте — отдельный класс задачи с тем же последовательным умножением, что SEQ. Для ускорения по потокам в этой
работе используются OMP и TBB.
