# Повышение контраста полутонового изображения посредством линейной растяжки гистограммы (ALL: MPI + STL)

**Вариант № 28**  

**Студент:** Клименко Владислав Сергеевич

**Группа:** 3823Б1ПР2  

**Преподаватель:** Сысоев Александр Владимирович, доцент

## 1. Введение

Алгоритм линейной растяжки гистограммы (contrast stretching) – один из базовых методов улучшения качества изображений.
Он преобразует входной диапазон яркостей пикселей [min_val, max_val] в полный диапазон [0, 255].

В работе представлена ALL-реализация – гибридный подход, комбинирующий MPI для межпроцессного взаимодействия и STL threads
для внутрипроцессного параллелизма.

Важное примечание: ALL-реализация не предназначена для прямого сравнения с OpenMP, TBB или STL на одном узле, т.к включает
дополнительные накладные расходы на MPI-коммуникации. Её основное преимущество раскрывается в распределённых средах.

## 2. Постановка задачи

Входные данные:
Вектор целых чисел input размера n, где каждый элемент – яркость пикселя (обычно от 0 до 255).

Выходные данные:
Вектор целых чисел output того же размера, содержащий значения после преобразования контрастности.

Преобразование:
Если max_val == min_val – изображение однородное, выход копируется с входа.
Иначе для каждого элемента:

        output[i] = ((input[i] - min_val) * 255) / (max_val - min_val)

Ограничения:

- n > 0

- Все значения лежат в диапазоне, допускающем целочисленное умножение без переполнения.

## 3. Базовый алгоритм (SEQ)

Последовательная версия выполняет два линейных прохода:

1. Поиск min_val и max_val с помощью std::ranges::minmax_element.

2. Преобразование каждого пикселя по формуле растяжки гистограммы.

Детальное описание приведено в seq/report.md.

## 4. ALL-алгоритм

4.1 Инициализация MPI и синхронизация размера.

        int rank = 0, size = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int total_size = 0;
        if (rank == 0) {
            total_size = static_cast<int>(input.size());
        }
        MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

4.2 Распределение данных (MPI_Scatterv).

        std::vector<int> recv_counts(size);
        std::vector<int> displs(size);

        int chunk = total_size / size;
        int remainder = total_size % size;
        for (int i = 0; i < size; i++) {
            recv_counts[i] = chunk + (i < remainder ? 1 : 0);
            displs[i] = (i == 0) ? 0 : displs[i - 1] + recv_counts[i - 1];
        }

        MPI_Scatterv(input.data(), recv_counts.data(), displs.data(), MPI_INT,
                    local_data.data(), local_size, MPI_INT, 0, MPI_COMM_WORLD);

4.3 Вспомогательная функция GetThreadRange.

        void GetThreadRange(size_t tid, size_t total, size_t num_t, size_t &begin, size_t &end) {
            size_t chunk = total / num_t;
            begin = tid * chunk;
            end = (tid == num_t - 1) ? total : begin + chunk;
        }

4.4 Параллельный поиск min/max (STL threads).

        std::pair<int, int> FindMinMaxSTL(const std::vector<int> &data, int num_threads) {
            // При пустых данных возвращаем предельные значения
            if (size == 0) return {INT_MAX, INT_MIN};
            // Запуск потоков для поиска локальных экстремумов
            for (int tid = 0; tid < num_threads; tid++) {
                threads[tid] = std::thread([&, tid]() {
                    // ... вычисление local_min, local_max
                });
            }
            // Объединение результатов
            return {global_min, global_max};
        }

4.5 Глобальная редукция (MPI_Allreduce).

        int global_min = 0, global_max = 0;
        MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

4.6 Преобразование контрастности (STL threads).

        void StretchContrastSTL(const std::vector<int> &input, std::vector<int> &output,
                                int min_val, int max_val, int num_threads) {
            // Параллельное преобразование каждого пикселя
            for (int tid = 0; tid < num_threads; tid++) {
                threads[tid] = std::thread([&, tid]() {
                    // ... применение формулы output[i] = ((input[i] - min) * 255) / (max - min)
                });
            }
        }

4.7 Сбор результатов (MPI_Gatherv).

        MPI_Gatherv(local_output.data(), local_size, MPI_INT,
                    recv_buffer.data(), recv_counts.data(), displs.data(), MPI_INT,
                    0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::ranges::copy(recv_buffer, output.begin());
        }

## 5. Implementation Details

- **Заголовочный файл:** `ops_all.hpp`.
- **Файл реализации:** `ops_all.cpp`.

## 6. Experimental Setup

- **Hardware/OS:**
  - CPU: Intel(R) Core(TM) i7-12650H
  - RAM: 16 ГБ
  - OS: Windows 10
- **Toolchain:**
  - IDE: Visual Studio Code
  - Компилятор: g++ 14.2.0
  - Тип сборки: Release

## 7. Результаты

### 7.1 Корректность

Функциональные тесты подтверждают правильность работы для:

- случайных массивов (50 + i % 101);

- однородных изображений (все пиксели равны);

- граничных размеров (1, 2, …, 1024).

Схема проверки:

- Для каждого тестового размера генерируется входной вектор по формуле 50 + (i % 101).

- Вычисляется эталонный результат по формуле преобразования.

- Результат последовательной реализации сравнивается с эталоном.

Результаты: все тесты пройдены. Расхождений между ожидаемым и фактическим выходом, а также между SEQ и OMP не обнаружено.

### 7.2 Производительность

Тесты на производительность успешно проходят проверку.

## 8. Выводы

ALL-реализация (MPI + STL) алгоритма повышения контрастности:

- полностью корректна – результаты совпадают с последовательной версией;

- обеспечивает масштабирование за пределы одного узла – основное преимущество перед чисто потоковыми решениями;

- демонстрирует гибкость – независимая настройка MPI и STL параллелизма;

- требует учёта накладных расходов на коммуникации (Scatterv/Gatherv/Allreduce).

Полученные результаты дают основу для дальнейшего расчёта ускорения и эффективности параллельных версий.

## 9. Список литературы

1. Лекции Сысоева Александра Владимировича по курсу «Параллельное программирование»

2. C++ Standard Library Reference – std::thread: <https://en.cppreference.com/w/cpp/thread/thread>

3. Документация по курсу: <https://learning-process.github.io/parallel_programming_course/ru/index.html>
