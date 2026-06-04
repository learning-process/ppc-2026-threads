# Умножение разреженных комплексных матриц (CCS) — SEQ

- Student: Иманов Сабутай Ширзад оглы, группа 3823Б1ПР5
- Technology: SEQ
- Variant: 7

## 1. Контекст

Данная работа посвящена **последовательной** реализации умножения двух разреженных комплексных матриц в формате CCS
(Compressed Column Storage). Последовательный вариант задаёт эталон корректности и базовое время $T_{seq}$ для сравнения
с OMP, TBB, STL и гибридной технологией ALL.

## 2. Постановка задачи

Краткая постановка: вычислить произведение $C = A \cdot B$ для матриц $A$, $B$ в CCS. Вход задачи — `InType =
std::tuple<CCS, CCS>`, выход — `OutType = CCS`. В `ValidationImpl` проверяется согласование размеров (`left.col_count ==
right.row_count`) и корректность структуры CCS (длины `col_start`, `row_index`, `nz`, допустимые индексы строк). При
слиянии ненулевых элементов отбрасываются значения с модулем не больше $10^{-14}$ (`kDropMagnitude` в
`all/src/ops_all.cpp`).

## 3. Базовый алгоритм

Алгоритм **столбцового** разреженного умножения (SpMM):

1. для каждого столбца $j$ правой матрицы $B$ собрать вклады из соответствующих столбцов $A$ (`BuildColumnFromRight`);
2. отсортировать пары (строка, значение) по номеру строки;
3. слить одинаковые строки и записать столбец в результирующую CCS (`CoalesceSortedPairs`).

Для SEQ все столбцы обрабатываются **в одном потоке** в функции `SpmmAbc`.

- **Время:** на каждый столбец $j$ — обход ненулевых $B$ и соответствующих столбцов $A$, сортировка буфера; в худшем
случае порядок растёт с числом ненулевых в столбце и размером буфера (сортировка $O(k \log k)$ для $k$ пар в буфере).
- **Память:** $O(\text{nnz}(C))$ для результата плюс буфер столбца (в коде `reserve(128)`).
- **Инварианты:** `col_start` монотонен; индексы строк в $[0, \text{row\_count})$; после слияния в столбце нет
дубликатов строк.

**Самая затратная часть:** `BuildColumnFromRight` и сортировка буфера для каждого столбца.

## 4. Схема распараллеливания

В режиме **SEQ** распараллеливание **не применяется**: вычисление выполняется одним потоком. Переменная окружения
`PPC_NUM_THREADS` на логику `SpmmAbc` не влияет; её изменение при замерах используется только для единообразия методики
курса (сравнение с OMP/TBB при том же числе потоков в окружении).

## 5. Детали реализации

Файлы: `seq/include/ops_seq.hpp`, `seq/src/ops_seq.cpp`, общее ядро — `all/src/ops_all.cpp`.

Каркас `Task` (как в примере курса):

| Метод | Роль в SEQ |
| ----- | ---------- |
| `ValidationImpl` | размеры $A$, $B$ и корректность CCS |
| `PreProcessingImpl` | без подготовки (`return true`) |
| `RunImpl` | вызов `BuildProductMatrix` → `SpmmAbc` |
| `PostProcessingImpl` | без постобработки (`return true`) |

Класс `SabutaySparseComplexCcsMultFixSEQ` делегирует умножение общему роутеру с типом `kSEQ`:

```45:63:tasks/sabutay_sparse_complex_ccs_mult_all/seq/src/ops_seq.cpp
void SabutaySparseComplexCcsMultFixSEQ::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out) {
  SabutaySparseComplexCcsMultAll::BuildProductMatrix(left, right, out, ppc::task::TypeOfTask::kSEQ);
}

bool SabutaySparseComplexCcsMultFixSEQ::RunImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  BuildProductMatrix(left, right, GetOutput());
  return true;
}
```

### Листинг: последовательный цикл по столбцам

```127:149:tasks/sabutay_sparse_complex_ccs_mult_all/all/src/ops_all.cpp
void SpmmAbc(const CCS &left, const CCS &right, CCS &out, double drop_magnitude) {
  out.row_count = left.row_count;
  out.col_count = right.col_count;
  out.col_start.assign(static_cast<std::size_t>(out.col_count) + 1U, 0);
  out.row_index.clear();
  out.nz.clear();
  if (out.col_count == 0) {
    return;
  }

  std::vector<std::pair<int, std::complex<double>>> buffer;
  buffer.reserve(128U);

  for (int j = 0; j < right.col_count; ++j) {
    BuildColumnFromRight(left, right, j, buffer);
    if (buffer.empty()) {
      out.col_start[static_cast<std::size_t>(j) + 1U] = static_cast<int>(out.nz.size());
      continue;
    }
    SortBufferByRow(buffer);
    CoalesceSortedPairs(buffer, out, drop_magnitude);
    out.col_start[static_cast<std::size_t>(j) + 1U] = static_cast<int>(out.nz.size());
  }
}
```

_Пояснение:_ внешний цикл по $j$ выполняется последовательно. Для каждого столбца буфер очищается, заполняется вкладами,
сортируется и сливается в глобальные массивы `row_index` и `nz` результата.

## 6. Проверка корректности

Корректность подтверждена функциональными тестами `tests/functional/main.cpp`: три случая (`case_id` 0, 1, 2). Для
каждого строится плотный эталон `Matmul`, результат CCS сравнивается с плотным видом (`CcsEqualsDenseView`), допуск
$10^{-12}$.

| `case_id` | Смысл теста |
| :-------: | ----------- |
| 0 | $3 \times 3$ и $3 \times 2$, вещественные значения |
| 1 | $2 \times 2$ и $2 \times 2$, комплексные фазы |
| 2 | $4 \times 1$ и $1 \times 3$, разреженная «цепочка» |

Тест `SabutaySparseComplexCcsMultFixSEQ` — в общем наборе `AddFuncTask` вместе с OMP, STL, TBB и ALL.

**Результаты func-тестов:** 3/3 **PASSED**.

| `case_id` | Суффикс GTest | Результат |
| :-------: | ------------- | :-------: |
| 0 | `..._seq_enabled_0` | PASSED |
| 1 | `..._seq_enabled_1` | PASSED |
| 2 | `..._seq_enabled_2` | PASSED |

Команда:

```powershell
$env:PPC_NUM_PROC = "1"
mpiexec -n 1 build\bin\ppc_func_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*seq_enabled*"
```

## 7. Экспериментальная среда

- **CPU:** Intel(R) Core(TM) Ultra 5 125H (14 ядер / 18 логических потоков).
- **RAM:** 32 ГБ.
- **ОС:** Windows 11.
- **Сборка:** Release (MSVC, Visual Studio 2022).
- **Размер задачи (perf):** случайные CCS `BuildRandomCcs(80, 90, 2027, 7)` и `BuildRandomCcs(90, 70, 4044, 6)` — см.
`tests/performance/main.cpp`, `SetUp()`.
- **Число прогонов внутри замера:** 5 (каркас `PerfAttr::num_running` в `modules/performance`).
- **Команда запуска (только SEQ, режим `task_run`):**

```powershell
$env:PPC_NUM_PROC = "1"
$env:PPC_NUM_THREADS = "4"   # для каждого из 1, 2, 4, 8
build\bin\ppc_perf_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*seq_enabled*"
```

## 8. Результаты

Каркас perf-теста запускает **`task_run`** (только `RunImpl`) и **`pipeline`** (полный цикл задачи).

**Базовое время для всего проекта:** $T_{seq} = 0.00001390$ с при $p = 1$. Это значение используется как знаменатель
ускорения $S$ в отчётах OMP, TBB, STL и ALL.

Ускорение ($S$) и эффективность ($E$) в таблице SEQ считаются относительно того же $T_{seq}$.

### `task_run`

| $p$ | Время (с) | $S$ | $E$ |
| :-: | :-------: | :-: | :-: |
| 1 | 0.00001390 | 1.00 | 100% |
| 2 | 0.00001416 | 0.98 | 49% |
| 4 | 0.00001538 | 0.90 | 23% |
| 8 | 0.00002128 | 0.65 | 8% |

_Комментарий:_ для последовательной реализации **отсутствие ускорения при росте $p$ ожидаемо**: вычислительное ядро
`SpmmAbc` не использует OpenMP/TBB. Небольшой разброс времени при разных $p$ объясняется фоновой нагрузкой ОС и
кэш-эффектами, а не параллельным ускорением алгоритма.

### Режим `pipeline`

Полный цикл задачи: `ValidationImpl` → `PreProcessingImpl` → `RunImpl` → `PostProcessingImpl`.

| $p$ | Время (с) |
| :-: | :-------: |
| 1 | 0.00002352 |
| 2 | 0.00008262 |
| 4 | 0.00002692 |
| 8 | 0.00004238 |

## 9. Выводы

Последовательная реализация задаёт **эталон** корректности и базовое время для сравнения технологий. Основная работа
сосредоточена в `SpmmAbc` (`all/src/ops_all.cpp`); класс SEQ — тонкая обёртка каркаса `Task`. Масштабирование по числу
потоков для SEQ не предусмотрено постановкой алгоритма; дальнейшее ускорение достигается в вариантах OMP и TBB за счёт
параллельной обработки столбцов.
