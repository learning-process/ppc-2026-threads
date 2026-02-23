# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников)

- Студент: Санников Иван Михайлович, Группа: 3823Б1ФИ2
- Технология: SEQ
- Вариант: 9

## 1. Введение
Численное вычисление многомерных интегралов является важной задачей вычислительной математики и широко применяется в статистике, физическом моделировании, машинном обучении и задачах оптимизации. Во многих прикладных ситуациях аналитическое вычисление интегралов невозможно или слишком трудоёмко, поэтому используются численные методы.

Одним из базовых подходов является метод прямоугольников, который строит приближение интеграла как сумму значений подынтегральной функции в узлах сетки, умноженную на объём элементарной ячейки разбиения. 

Цель данной лабораторной работы — реализовать вычисление многомерных интегралов с использованием метода прямоугольников в последовательном варианте SEQ и провести проверку корректности результата на наборе тестовых задач.

## 2. Постановка задачи

Входные данные: `InType` — кортеж вида:

- `std::function<double(const std::vector<double>&)> func` — подынтегральная функция `f(x)`.
- `std::vector<std::pair<double, double>> borders` — границы интегрирования по каждому измерению.
- `int n` — число разбиений на каждом измерении.

Выходные данные: значение приближённого интеграла типа `double`.

Для области интегрирования \[a, b\] вычислить приближение интеграла методом прямоугольников по равномерной
сетке.

## 3. Базовый алгоритм (seq)

Последовательный алгоритм:

- Определяется размерность.
- Для каждого измерения вычисляется шаг: `h_i = (b_i - a_i) / n`.
- Вычисляется объём одной ячейки: `V` равен произведению всех `h_i`.
- Выполняется перебор всех ячеек сетки `(k_1, ..., k_d)`, где `k_i` принадлежит \[0, n - 1\].
- Для каждой ячейки вычисляется центр: `x_i = a_i + (k_i + 0.5) * h_i`.
- Суммируется `f(x)` по всем ячейкам и умножается на `V`.

```cpp
bool SannikovIIntegralsRectangleMethodSEQ::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const std::size_t dim = borders.size();
  std::vector<double> h(dim);
  double cell_v = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    const double left_border = borders[i].first;
    const double right_border = borders[i].second;
    h[i] = (right_border - left_border) / static_cast<double>(n);
    if (!(h[i] > 0.0) || !std::isfinite(h[i])) {
      return false;
    }
    cell_v *= h[i];
  }
  std::vector<int> idx(dim, 0);
  std::vector<double> x(dim);
  double sum = 0.0;
  while (true) {
    for (std::size_t i = 0; i < dim; ++i) {
      const double left_border = borders[i].first;
      x[i] = left_border + ((static_cast<double>(idx[i]) + 0.5) * h[i]);
    }
    const double fx = func(x);
    if (!std::isfinite(fx)) {
      return false;
    }
    sum += fx;

    if (!NextIndex(idx, dim, n)) {
      break;
    }
  }
  GetOutput() = sum * cell_v;
  return std::isfinite(GetOutput());
}
```

## 4. Описание параллельного алгоритма 
Данная работа не подразумевает использования и реализацию параллельного алгоритма.

## 5. Experimental Setup

- Hardware/OS: Intel i9 13900KF, 24 ядра, RAM: 32Gb, OS: Windows 11
- Toolchain: Cmake 3.28.3, g++ (Ubuntu 14.2.0 x86_64), Docker-контейнер, Режим сборки: Release.
- Data: Для тестов на производительность использовалась равномерная сетка разбиения области интегрирования
  \[0, π\] × \[0, π\] с числом разбиений 10000 по каждому измерению, где в каждой точке вычислялась функция
  `f(x, y) = sin(x) sin(y)` методом прямоугольников.

## 6. Результаты

### 6.1 Корректность 

Корректность работы алгоритма проверена с использованием Google Test. В тестах использовались интеграл от
константы, линейные и квадратичные функции, произведения переменных, тригонометрические и экспоненциальные
функции.

### 6.2 Производительность


| Mode        |  Size of the split | Time, s |
|-------------|--------------------|---------|
| seq         | 500                | 0.0022  |
| seq         | 1000               | 0.0088  |
| seq         | 5000               | 0.2290  |
| seq         | 10000              | 0.9364  |

Так как алгоритм выполняет `n^dim` вычислений, наблюдается квадратичная зависимость времени работы от числа
разбиений (при увеличении числа разбиений в 2 раза время увеличивается примерно в 4 раза).
## 7. Выводы

В рамках лабораторной работы реализован последовательный алгоритм вычисления многомерного интеграла методом
прямоугольников. Программа поддерживает произвольную размерность области интегрирования, заданную через список
границ, и использует равномерное разбиение по каждому измерению. Корректность результата подтверждена
функциональными тестами на наборе функций с известными аналитическими интегралами.

## 8. Литература
1. Керченский государственный морской технологический университет.  Методы прямоугольников - https://studfile.net/preview/10730514/page:2/
2. Parallel Programming 2025-2026 - https://disk.yandex.com/d/KVEpMAlUWXwg8Q

## 9. Приложение

```cpp

namespace {
bool NextIndex(std::vector<int>& idx, std::size_t dim, int n) {
  for (std::size_t pos = 0; pos < dim; ++pos) {
    ++idx[pos];
    if (idx[pos] < n) {
      return true;
    }
    idx[pos] = 0;
  }
  return false;
}

}  // namespace

bool SannikovIIntegralsRectangleMethodSEQ::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const std::size_t dim = borders.size();

  std::vector<double> h(dim);
  double cell_v = 1.0;

  for (std::size_t i = 0; i < dim; ++i) {
    const double left_border = borders[i].first;
    const double right_border = borders[i].second;

    h[i] = (right_border - left_border) / static_cast<double>(n);
    if (!(h[i] > 0.0) || !std::isfinite(h[i])) {
      return false;
    }

    cell_v *= h[i];
  }

  std::vector<int> idx(dim, 0);
  std::vector<double> x(dim);

  double sum = 0.0;

  while (true) {
    for (std::size_t i = 0; i < dim; ++i) {
      const double left_border = borders[i].first;
      x[i] = left_border + ((static_cast<double>(idx[i]) + 0.5) * h[i]);
    }

    const double fx = func(x);
    if (!std::isfinite(fx)) {
      return false;
    }

    sum += fx;

    if (!NextIndex(idx, dim, n)) {
      break;
    }
  }

  GetOutput() = sum * cell_v;

  return std::isfinite(GetOutput());
}

```