# Многомерное интегрирование Монте-Карло — SEQ

- **Student:** Шилин Никита Дмитриевич, группа 3823Б1ПР1
- **Technology:** SEQ
- **Variant:** 12

---

## 1. Контекст

Задача — оценить значение определённого интеграла по осевому параллелепипеду
\([\mathrm{lower}, \mathrm{upper}]\) методом Монте-Карло с **квази-случайной**
последовательностью Кронекера (детерминированные иррациональные сдвиги
\(\alpha_d\)). Последовательная реализация задаёт **эталон корректности и эталон
времени** для всех остальных backend-ов; именно её время используется как
знаменатель в формулах ускорения \(S = T_{\mathrm{seq}}/T_p\) в потоковых
(`omp`/`tbb`/`stl`) и гибридной (`all`) версиях.

## 2. Постановка задачи

**Вход.** Кортеж `InType = std::tuple<std::vector<double>, std::vector<double>, int, FuncType>`,
определённый в `common/include/common.hpp`, — это `(lower, upper, n, func_type)`:
границы по осям, число выборок \(n\) и тип подынтегральной функции из
[`IntegrandFunction`](../common/include/common.hpp).

**Выход.** `OutType = double` — оценка интеграла \(I \approx V \cdot \tfrac{1}{n}\sum_i f(x^{(i)})\),
где \(V = \prod_d (\mathrm{upper}_d - \mathrm{lower}_d)\).

**Ограничения.** Они проверяются ровно в `ValidationImpl()`:

```16:34:tasks/shilin_n_monte_carlo_integration/seq/src/ops_seq.cpp
bool ShilinNMonteCarloIntegrationSEQ::ValidationImpl() {
  const auto &[lower, upper, n, func_type] = GetInput();
  if (lower.size() != upper.size() || lower.empty()) {
    return false;
  }
  if (n <= 0) {
    return false;
  }
  for (size_t i = 0; i < lower.size(); ++i) {
    if (lower[i] >= upper[i]) {
      return false;
    }
  }
  if (func_type < FuncType::kConstant || func_type > FuncType::kSinProduct) {
    return false;
  }
  constexpr size_t kMaxDimensions = 10;
  return lower.size() <= kMaxDimensions;
}
```

Граничные случаи (пустой вектор, `n ≤ 0`, перевёрнутые границы, недопустимый
`FuncType`, размерность > 10) приводят к `false` и отказу инфраструктуры курса
запускать `RunImpl`.

**Корректность.** Эталоном служит аналитический интеграл
[`IntegrandFunction::AnalyticalIntegral`](../common/include/common.hpp);
тест считает прогон корректным, если
\(|\hat I - I_{\mathrm{exact}}| \le \max\!\bigl(10\,V/\sqrt{n},\,10^{-2}\bigr)\) —
это согласовано со стандартной оценкой ошибки Монте-Карло \(O(1/\sqrt{n})\).

## 3. Базовый алгоритм

Используется детерминированная **последовательность Кронекера**:
для каждой оси \(d\) фиксируется иррациональная константа
\(\alpha_d = \{\sqrt{p_d}\}\) — дробная часть корня из простого числа
(\(p_d \in \{2,3,5,7,11,13,17,19,23,29\}\)). Координата вычисляется как

\[
t_d^{(i)} = \{0.5 + (i+1)\,\alpha_d\}, \qquad x_d^{(i)} = \mathrm{lower}_d + (\mathrm{upper}_d-\mathrm{lower}_d)\,t_d^{(i)}.
\]

В `RunImpl` накопление выполняется итерационно: в `current[d]` хранится
текущая дробная часть, которая обновляется добавлением \(\alpha_d\) и
последующим вычитанием единицы при переполнении.

```46:84:tasks/shilin_n_monte_carlo_integration/seq/src/ops_seq.cpp
bool ShilinNMonteCarloIntegrationSEQ::RunImpl() {
  auto dimensions = static_cast<int>(lower_bounds_.size());

  // Quasi-random Kronecker sequence: alpha_d = fractional part of sqrt(prime_d)
  const std::vector<double> alpha = { /* sqrt(2), sqrt(3), ... */ };

  std::vector<double> current(dimensions, 0.5);
  std::vector<double> point(dimensions);
  double sum = 0.0;

  for (int i = 0; i < num_points_; ++i) {
    for (int di = 0; di < dimensions; ++di) {
      current[di] += alpha[di];
      if (current[di] >= 1.0) {
        current[di] -= 1.0;
      }
      point[di] = lower_bounds_[di] + ((upper_bounds_[di] - lower_bounds_[di]) * current[di]);
    }
    sum += IntegrandFunction::Evaluate(func_type_, point);
  }
  ...
  GetOutput() = volume * sum / static_cast<double>(num_points_);
```

**Сложность.** \(O(n \cdot d)\) по времени и \(O(d)\) по памяти (помимо
`InType/OutType`); для perf-теста курса \(n = 10^7\), \(d = 3\) — это
\(\approx 3 \cdot 10^7\) операций с плавающей точкой плюс \(10^7\) вызовов
`IntegrandFunction::Evaluate`.

## 4. Детали реализации

- Файлы: [`seq/include/ops_seq.hpp`](include/ops_seq.hpp), [`seq/src/ops_seq.cpp`](src/ops_seq.cpp).
- Класс `ShilinNMonteCarloIntegrationSEQ` — наследник `BaseTask = ppc::task::Task<InType,OutType>`,
  возвращает `TypeOfTask::kSEQ` и переопределяет четыре метода каркаса
  `Validation/PreProcessing/Run/PostProcessing`.
- В `PreProcessingImpl` входной кортеж распаковывается в поля класса
  (`lower_bounds_`, `upper_bounds_`, `num_points_`, `func_type_`); это разделяет
  валидацию и собственно вычисление и снимает повторное structured-binding
  внутри `RunImpl`.
- `PostProcessingImpl` — `return true`: дополнительная обработка не нужна,
  потому что `RunImpl` пишет результат сразу в `GetOutput()`.

## 5. Проверка корректности

Общий тестовый набор курса определён в
[`tests/functional/main.cpp`](../tests/functional/main.cpp) — 8 параметризованных
кейсов, прогоняемых для всех backend-ов через `ppc::util::AddFuncTask<...>`:

| Кейс | \(d\) | \(n\) | Функция |
| --- | ---: | ---: | --- |
| `kLinear_1D` | 1 | 10 000 | `kLinear` |
| `kSumSquares_1D` | 1 | 10 000 | `kSumSquares` |
| `kConstant_2D` | 2 | 50 000 | `kConstant` |
| `kLinear_2D` | 2 | 50 000 | `kLinear` |
| `kProduct_2D` | 2 | 50 000 | `kProduct` |
| `kSinProduct_2D` | 2 | 50 000 | `kSinProduct` |
| `kLinear_3D` | 3 | 100 000 | `kLinear` |
| `kProduct_3D` | 3 | 100 000 | `kProduct` |

Эти же кейсы используют OMP/TBB/STL/ALL — это и есть реализация
требования преподавателя «функциональные тесты унифицированы».
Локальный прогон (Apple M4 Max, build-local) — все 8 кейсов SEQ проходят.

## 6. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | Apple M4 Max, 16 ядер (12P + 4E) |
| RAM | 64 GiB |
| OS | macOS 26.3.1 (build 25D771280a) |
| Компилятор | Apple clang 17.0.0 (`/usr/bin/c++`) |
| Сборка | `Release`, `Unix Makefiles`, `build-local` |
| Каркас perf | `ppc::util::BaseRunPerfTests<InType, OutType>` |
| Размер задачи | \(n = 10^7\), \([0,1]^3\), `FuncType::kSumSquares` |
| Повторов | 5 запусков, отчётно — медиана |

Команда сборки и запуска:

```bash
cmake -S . -B build-local -G "Unix Makefiles" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/c++ \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
  -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON
cmake --build build-local -j 16
./build-local/bin/ppc_perf_tests --gtest_filter='*shilin_n_monte_carlo_integration_seq*'
```

## 7. Результаты

| Режим | Время, с (медиана 5) | \(S\) vs SEQ | Комментарий |
| --- | ---: | ---: | --- |
| `task_run` | **0,030436** | 1,00 | базовое значение для расчёта \(S\) и \(E\) |
| `pipeline` | **0,030642** | 1,00 | `task_run` и `pipeline` для SEQ совпадают в пределах шума (~0,7%): для последовательной версии конвейер каркаса и сам `RunImpl` доминируют примерно одинаково |

Самая дорогая часть `RunImpl` — внутренний цикл `for (int i = 0; i < num_points_; ++i)` с
вызовом `IntegrandFunction::Evaluate`; на `kSumSquares_3D` это `O(n·d)`
плавающих умножений и сложений плюс вычисление `floor` для дробной части,
которое полностью векторизуется компилятором.

## 8. Выводы

- SEQ-версия используется как **знаменатель** для всех остальных моделей
  параллелизма; время измеряется в одинаковом каркасе perf-тестов курса.
- Алгоритм по структуре — **embarrassingly parallel** по индексу \(i\):
  единственный data dependency — глобальная сумма, которая редуцируется в
  одну переменную (это хорошо ложится на OMP `reduction`, TBB `parallel_reduce`,
  ручную редукцию в STL и `MPI_Allreduce` в ALL).
- Стабильное `task_run`-время **≈0,0304 с** на M4 Max — отправная точка
  для оценки реального ускорения остальных backend-ов в их соответствующих
  отчётах.

## 9. Источники

- Документация курса PPC, репозиторий `ppc-2026-threads`.
- Niederreiter H. *Random Number Generation and Quasi-Monte Carlo Methods.*
- C++ reference: [`std::tuple`](https://en.cppreference.com/w/cpp/utility/tuple).

## 10. Чек-лист

- [x] Постановка задачи и ограничения зафиксированы.
- [x] Базовый алгоритм описан, сложность указана (`O(n·d)`).
- [x] Корректность подтверждена 8 функциональными кейсами и сравнением с
      аналитическим интегралом.
- [x] Среда замеров и команды воспроизведения указаны.
- [x] Базовое `task_run`-время есть (медиана 5 прогонов).
- [x] Этот файл — `seq/report.md`, лежит ровно в `tasks/<task>/seq/`.
