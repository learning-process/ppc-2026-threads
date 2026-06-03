# Умножение разреженных матриц комплексного типа в формате CRS — SEQ

- Student: Климович Виктор Олегович, group 3823Б1ПР4
- Technology: SEQ
- Variant: 6

## 1. Контекст

Последовательная реализация служит эталоном корректности и базой для расчёта
ускорений всех остальных backend-ов (OMP, TBB, STL, ALL). От её точности
зависит методическая обоснованность сравнения технологий: если в SEQ есть
скрытая ошибка или нестабильная численная агрегация, любое полученное
«ускорение» будет сравнением неверного результата с неверным результатом.
Поэтому здесь подробно фиксируется алгоритм, инварианты и крайние случаи.

## 2. Постановка задачи

На вход поступают две разреженные матрицы `A` и `B` с элементами типа
`std::complex<double>`, хранимые в формате Compressed Row Storage (CRS):
для каждой матрицы заданы вектора `row_offsets`, `col_indices`, `data` и поля
`n_rows`, `n_cols` (см. [common/include/common.hpp](../common/include/common.hpp)).
Требуется вычислить разреженное произведение `C = A * B` в том же формате CRS.

Ограничения на вход:

- `A.n_cols == B.n_rows` (проверяется в `ValidationImpl`).
- `row_offsets[0] == 0`, `row_offsets[i] <= row_offsets[i+1]`, длина массива
  `row_offsets` равна `n_rows + 1`.
- Элементы внутри строки CRS отсортированы по возрастанию `col_indices`.

Дополнительный численный критерий: значения результата, у которых
`|Re| <= kZeroDropTol = 1e-12` И `|Im| <= kZeroDropTol`, отбрасываются — это
обычная санитизация для накопления почти-нулей в разреженных алгоритмах.

## 3. Базовый алгоритм

Используется row-wise алгоритм Густавсона со sparse accumulator (SPA). Для
каждой строки `i` матрицы `A` обходятся её ненулевые элементы `(i, k)`, и для
каждого из них — ненулевые элементы строки `k` матрицы `B` `(k, j)`; вклад
`A[i,k] * B[k,j]` накапливается в плотном буфере `spa[j]` размера `B.n_cols`.

Чтобы не сбрасывать весь `spa` после каждой строки, используется sentinel-трюк:
вектор `touched_by_row[j]` хранит индекс «последней строки, в которой `spa[j]`
был задействован». Условие «первая запись в `spa[j]` для текущей строки `i`»
проверяется как `touched_by_row[j] != i`. Это даёт амортизированную сложность
**O(nnz(A) + Σ_i Σ_{k ∈ row(A,i)} nnz(B[k, :]))** во времени и
**O(B.n_cols)** в дополнительной памяти на SPA + `touched_by_row` + `touched_cols`.

Инварианты:

- После обработки строки `i` все `spa[j]` для затронутых столбцов обнулены,
  а `touched_by_row[j]` хранит значение `i` (он не сбрасывается — нужен лишь
  для будущих сравнений; начальное значение `-1` гарантирует, что для `i = 0`
  ни одна ячейка не считается «затронутой»).
- `touched_cols` хранит только реально записанные столбцы текущей строки;
  после сортировки порядок CRS соблюдён.
- `result.row_offsets[i+1] = result.row_offsets[i] + kept`, где `kept` — число
  оставленных после фильтрации почти-нулей значений.

## 4. Детали реализации

Файлы: [seq/include/ops_seq.hpp](include/ops_seq.hpp),
[seq/src/ops_seq.cpp](src/ops_seq.cpp).

- `ValidationImpl`: проверяет согласованность размерностей
  `A.n_cols == B.n_rows`.
- `PreProcessingImpl`: пустой — данные уже в нужном формате, выделение
  результата делается внутри `RunImpl`.
- `RunImpl`: вызывает `MultiplyCrs(lhs, rhs)`.
- `PostProcessingImpl`: пустой — формат CRS строится напрямую во время
  основного цикла.

Ключевые фрагменты — две вспомогательные функции в анонимном namespace:

```cpp
// File: seq/src/ops_seq.cpp
void AccumulateRow(const CrsMatrix &lhs, const CrsMatrix &rhs, int row,
                   std::vector<Cplx> &spa,
                   std::vector<int> &touched_by_row,
                   std::vector<int> &touched_cols) {
  touched_cols.clear();
  for (int lp = lhs.row_offsets[row]; lp < lhs.row_offsets[row + 1]; ++lp) {
    const int k = lhs.col_indices[lp];
    const Cplx a_ik = lhs.data[lp];
    for (int rq = rhs.row_offsets[k]; rq < rhs.row_offsets[k + 1]; ++rq) {
      const int j = rhs.col_indices[rq];
      if (touched_by_row[j] != row) {
        touched_by_row[j] = row;
        touched_cols.push_back(j);
        spa[j] = a_ik * rhs.data[rq];
      } else {
        spa[j] += a_ik * rhs.data[rq];
      }
    }
  }
}
```

В этой функции `touched_by_row[j] != row` — это и есть sentinel-трюк: для
текущей строки `row` любая «свежая» запись в `spa[j]` помечается флагом, и
при последующих обращениях к тому же `j` мы знаем, что значение можно
накапливать.

```cpp
void FinalizeRow(CrsMatrix &result, std::vector<Cplx> &spa,
                 std::vector<int> &touched_cols, int row) {
  std::ranges::sort(touched_cols);
  int kept = 0;
  for (const int j : touched_cols) {
    const Cplx v = spa[j];
    spa[j] = Cplx(0.0);
    if (std::abs(v.real()) > kZeroDropTol || std::abs(v.imag()) > kZeroDropTol) {
      result.col_indices.push_back(j);
      result.data.push_back(v);
      ++kept;
    }
  }
  result.row_offsets[row + 1] = result.row_offsets[row] + kept;
}
```

Здесь столбцы сортируются (для соблюдения инварианта CRS), `spa[j]`
обнуляется сразу после чтения (это и есть «ленивая» очистка только затронутых
ячеек), а малые по модулю значения отфильтровываются по `kZeroDropTol`.

## 5. Проверка корректности

Используется functional-тест
[tests/functional/main.cpp](../tests/functional/main.cpp): для пар размеров
`(rows_a, cols_a_rows_b, cols_b)` ∈ `kTestParams` (10 наборов, от `1×1×1` до
`11×3×4` и `7×12×9`) генерируется псевдослучайная плотная пара матриц `A`,
`B`, считается эталон `DenseMultiply(A, B)`, а результат SEQ переводится
обратно в плотный вид через `ToDense` и сравнивается поэлементно с допуском
`kCompareTol = 1e-9`. Генератор `GenerateDense` даёт ≈25% ненулей, что
обеспечивает покрытие как «вертикально-узких», так и «горизонтально-узких»
случаев.

Все 10 наборов проходят на SEQ-реализации без расхождений.

## 6. Экспериментальная среда

- CPU: Intel Core i7-11800H @ 2.30 GHz (8 физических ядер / 16 аппаратных потоков, Tiger Lake-H)
- RAM: 16 GB DDR4-3200
- OS: Windows 11 Pro 22H2 (build 22631)
- Compiler: MSVC 19.44.35211 (Visual Studio 2022)
- Build type: Release
- Размеры задач: 10 наборов functional, `1500×1500` (≈2 ненуля на строку) для
  performance.

Команды:

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --parallel

set PPC_NUM_THREADS=1
scripts/run_tests.py --running-type=threads --counts 1
scripts/run_tests.py --running-type=performance
```

## 7. Результаты

Baseline для `p = 1` на матрице `1500×1500`:

| size      | mode     | median time, s |
| --------- | -------- | -------------: |
| 1500×1500 | task     |          0.152 |
| 1500×1500 | pipeline |          1.498 |

Эти значения служат знаменателями `T_seq` для расчёта ускорения
`S(p) = T_seq / T(p)` во всех остальных локальных отчётах.

Самым дорогим фрагментом ожидаемо является вложенный цикл обхода
`row(A) × row(B)` со случайным доступом в `spa` по `j`; именно он определяет
кэш-локальность и будет лимитировать масштабируемость параллельных версий.

## 8. Выводы

Реализация SEQ — однопоточный Густавсон с одним общим SPA, который
переиспользуется между строками за счёт sentinel-обнуления. Сложность
оптимальна для row-wise сценария, корректность подтверждена сравнением с
плотным эталоном на всех функциональных тестах. Эта версия принята как
`T_seq = 0.152 s` (task) и `T_seq = 1.498 s` (pipeline) для всех последующих
сравнений.
