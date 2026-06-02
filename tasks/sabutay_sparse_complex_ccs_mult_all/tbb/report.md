# Умножение разреженных комплексных матриц (CCS) — TBB

- Student: Иманов Сабутай Ширзад оглы, группа 3823Б1ПР5
- Technology: TBB
- Variant: 7

## 1. Контекст

В этой версии для параллельной работы используется библиотека **oneTBB** (`parallel_for`). Столбцы результата считаются
параллельно, как в OMP. Эталон SEQ: $T_{seq} = 0.00001390$ с.

## 2. Постановка задачи

Нужно вычислить произведение $C = A \cdot B$ для двух разреженных **комплексных** матриц в формате **CCS**.

- **Вход:** `InType = std::tuple<CCS, CCS>`.
- **Выход:** `OutType = CCS`.
- **Проверка:** `ValidationImpl` класса `SabutaySparseComplexCcsMultFixTBB` — размеры матриц и корректность CCS.
- **Порог слияния:** $10^{-14}$ по модулю элемента.
- **Сравнение в тестах:** с плотным умножением, допуск $10^{-12}$.
- **Эталон времени:** $T_{seq} = 0.00001390$ с.

## 3. Базовый алгоритм

Столбцовое умножение CCS (как в SEQ):

1. `BuildColumnFromRight` — собрать столбец результата для заданного $j$;
2. `SortBufferByRow` — сортировка по строке;
3. `CoalesceSortedPairs` — слияние и запись в CCS.

В **TBB** цикл по $j$ выполняется через `oneapi::tbb::parallel_for` в `BuildProductMatrixTbb`. Каждый столбец пишется в
свой локальный вектор, затем столбцы склеиваются в одну матрицу **последовательно**.

## 4. Схема распараллеливания

- `tbb::parallel_for` делит диапазон столбцов $[0, \text{col\_count})$ на части.
- Для каждого столбца результат пишется в отдельные векторы `local_row_index[j]` и `local_nz[j]` — гонок данных нет.
- После `parallel_for` один поток последовательно склеивает столбцы в итоговую матрицу.

Число потоков — через `PPC_NUM_THREADS`.

**Параметры TBB:**

| Элемент | В проекте |
| ------- | --------- |
| Примитив | `oneapi::tbb::parallel_for` |
| Диапазон | `blocked_range<int>(0, right.col_count)` |
| Размер зерна (grainsize) | по умолчанию (авторазбиение runtime) |
| Partitioner | по умолчанию |
| Гонки | нет: столбец $j$ пишет только в `local_row_index[j]`, `local_nz[j]` |

## 5. Детали реализации

Файлы: `tbb/src/ops_tbb.cpp`, функция `BuildProductMatrixTbb` в `all/src/ops_all.cpp`.

```45:63:tasks/sabutay_sparse_complex_ccs_mult_all/tbb/src/ops_tbb.cpp
void SabutaySparseComplexCcsMultFixTBB::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out) {
  SabutaySparseComplexCcsMultAll::BuildProductMatrix(left, right, out, ppc::task::TypeOfTask::kTBB);
}

bool SabutaySparseComplexCcsMultFixTBB::RunImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  BuildProductMatrix(left, right, GetOutput());
  return true;
}
```

### Листинг: `parallel_for` по столбцам

```213:231:tasks/sabutay_sparse_complex_ccs_mult_all/all/src/ops_all.cpp
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, right.col_count),
                            [&](const oneapi::tbb::blocked_range<int> &range) {
    std::vector<std::pair<int, std::complex<double>>> buffer;
    buffer.reserve(128U);
    for (int jcol = range.begin(); jcol < range.end(); ++jcol) {
      BuildColumnFromRight(left, right, jcol, buffer);
      if (!buffer.empty()) {
        SortBufferByRow(buffer);
        CCS tmp;
        CoalesceSortedPairs(buffer, tmp, kDropMagnitude);
        local_row_index[static_cast<std::size_t>(jcol)] = std::move(tmp.row_index);
        local_nz[static_cast<std::size_t>(jcol)] = std::move(tmp.nz);
      }
      ...
    }
  });
```

_Пояснение:_ TBB сам распределяет столбцы по потокам. Каждый столбец хранится отдельно, потом всё объединяется в одну
CCS.

## 6. Проверка корректности

`tests/functional/main.cpp`, класс `SabutaySparseComplexCcsMultFixTBB` — сравнение с плотным умножением.

**Результаты func-тестов:** 3/3 **PASSED**.

| `case_id` | Суффикс GTest | Результат |
| :-------: | ------------- | :-------: |
| 0 | `..._tbb_enabled_0` | PASSED |
| 1 | `..._tbb_enabled_1` | PASSED |
| 2 | `..._tbb_enabled_2` | PASSED |

```powershell
$env:PPC_NUM_PROC = "1"
mpiexec -n 1 build\bin\ppc_func_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*tbb_enabled*"
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
build\bin\ppc_perf_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*tbb_enabled*"
```

## 8. Результаты

Каркас perf запускает **`task_run`** и **`pipeline`**. $T_{seq} = 0.00001390$ с.

### `task_run`

| $p$ | Время (с) | $S$ | $E$ |
| :-: | :-------: | :-: | :-: |
| 1 | 0.00008308 | 0.17 | 17% |
| 2 | 0.00011736 | 0.12 | 6% |
| 4 | 0.00011596 | 0.12 | 3% |
| 8 | 0.00014240 | 0.10 | 1% |

_Комментарий:_ при одном потоке TBB уже медленнее SEQ — библиотека тратит время на организацию параллельной работы. На
маленькой perf-задаче ускорения нет.

### `pipeline`

| $p$ | Время (с) |
| :-: | :-------: |
| 1 | 0.00021300 |
| 2 | 0.00024880 |
| 4 | 0.00023686 |
| 8 | 0.00027230 |

## 9. Выводы

TBB удобен, когда много независимых частей работы (здесь — столбцы). Для матриц из perf-теста выгоднее простой SEQ.
Результаты совпадают с эталоном.
