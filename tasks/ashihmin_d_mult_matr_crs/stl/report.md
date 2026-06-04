# Умножение разреженных матриц в формате CRS - STL

* Student: Ашихмин Д., group 3823Б1ФИ2
* Technology: STL
* Variant: 4

## 1. Introduction

STL-версия использует стандартные механизмы многопоточности C++, такие как
`std::async` и `std::future`. Эта реализация демонстрирует, как задачу умножения
разреженных матриц можно эффективно распараллелить, используя исключительно
средства стандартной библиотеки (ISO C++), что обеспечивает максимальную переносимость
кода без зависимости от OpenMP или TBB.

## 2. Problem Statement

Задача совпадает с SEQ-версией: реализовать умножение двух разреженных матриц A и B,
представленных в строковом формате (CRS).
Вход: InType = std::pair<CRSMatrix, CRSMatrix>.
Выход: OutType = CRSMatrix, результат C = A * B.
Ограничения: matrix_a.cols == matrix_b.rows, элементы типа double.

## 3. Baseline Algorithm (Sequential)

Baseline описан в seq/report.md. Последовательный алгоритм выполняет обход строк матрицы A,
вычисляет ненулевые значения результирующей строки через вспомогательный аккумулятор и сохраняет
их в структуру CRS.

## 4. Parallelization Scheme

Для распределения нагрузки используется декомпозиция по строкам матрицы A. Весь диапазон
строк делится на блоки (chunks), количество которых соответствует числу аппаратных потоков:

* **Определение потоков:** thread_count берется из `std::thread::hardware_concurrency()`.
* **Разбиение на блоки:** Каждому потоку назначается диапазон строк `[start_row, end_row)`.
* **Запуск задач:** Для каждого блока вызывается `std::async` с политикой `std::launch::async`,
что гарантирует выполнение в отдельном потоке.
* **Вычисления:** Внутри каждой задачи выполняется цикл по назначенным строкам.
Для каждой строки вызывается метод `MultiplyRow`, использующий локальный `std::map` для накопления значений.
* **Синхронизация:** Главный поток сохраняет объекты `std::future` в вектор и дожидается
завершения всех задач через вызов `.get()`.
* **Безопасность:** Гонок данных нет, так как каждый поток пишет результаты только в свои
индексы векторов `local_cols` и `local_vals`.

## 5. Implementation Details

* Файлы: stl/include/ops_stl.hpp, stl/src/ops_stl.cpp.
* Класс: AshihminDMultMatrCrsSTL.
* **MultiplyRow:** Логика вычисления одной строки вынесена в статический метод для снижения
когнитивной сложности и удовлетворения требований Clang-Tidy.
* **Память:** Дополнительная память используется для хранения промежуточных локальных векторов
каждой строки перед финальной сборкой.

## 6. Experimental Setup

Аппаратное обеспечение:

* CPU: AMD Ryzen 5 3500X (3.60 GHz, 6 ядер / 6 потоков)
* RAM: 16 ГБ
* OS: Windows 11 / Linux (CI)

Инструменты:

* Сборка: CMake
* Компилятор: GCC / Clang / MSVC (стандарт C++20)

## 7. Results and Discussion

### 7.1 Correctness

Корректность проверялась функциональными тестами из tests/functional/main.cpp.
Вычисленные значения сравниваются с результатами последовательной версии для различных
типов матриц (единичные, прямоугольные, случайные) с точностью 1e-10.

### 7.2 Performance

Mode | Count | Time, s | Speedup | Efficiency
--- | --- | --- | --- | ---
seq | 1 | 0.852412 | 1.00 | N/A
stl | 2 | 0.463267 | 1.84 | 92.00%
stl | 4 | 0.247794 | 3.44 | 86.00%
stl | 6 | 0.181245 | 4.70 | 78.33%

## 8. Conclusions

Реализация на базе STL (`std::async`) показала производительность, сопоставимую с
OpenMP. Использование chunks (блоков строк) вместо создания потока на каждую строку
позволило минимизировать накладные расходы на управление задачами.

## 9. References

1. ISO C++ Standard Library Documentation: std::thread, std::async.
2. Anthony Williams. C++ Concurrency in Action.
3. OpenMP Architecture Review Board. OpenMP API.
4. Microsoft MPI Documentation.

## Appendix (Optional)

Фрагмент RunImpl() с использованием std::async:

```cpp
for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
  int start_row = thread_idx * chunk_size;
  int end_row = std::min(start_row + chunk_size, rows_a);
  if (start_row >= end_row) break;

  futures.push_back(std::async(std::launch::async, [=, &matrix_a, &matrix_b, &local_cols, &local_vals] {
    for (int i = start_row; i < end_row; ++i) {
      MultiplyRow(i, matrix_a, matrix_b, local_cols[i], local_vals[i]);
    }
  }));
}
for (auto &fut : futures) fut.get();
