# Поразрядная сортировка для целых чисел с простым слиянием

- **Студент:** Вотинцев Д. С.
- **Группа:** 3823Б1ФИ3
- **Технология:** ALL (MPI + OpenMP)
- **Вариант:** 17

## 1. Введение

Гибридные модели программирования, сочетающие MPI (Message Passing Interface) и OpenMP, позволяют эффективно
использовать вычислительные ресурсы современных кластерных систем. В такой схеме MPI отвечает за распределение
данных между узлами (процессами), а OpenMP — за параллелизм внутри каждого узла на уровне ядер процессора.
Поразрядная сортировка (Radix Sort), обладая линейной сложностью, хорошо поддается такому двухуровневому
распараллеливанию.

## 2. Постановка задачи

Реализовать гибридный алгоритм поразрядной сортировки для 32-битных целых чисел.

- **Входные данные**: Вектор `std::vector<int32_t>`.
- **Выходные данные**: Отсортированный вектор на нулевом процессе.
- **Особенности**: Использование MPI для передачи сообщений и OpenMP для ускорения локальных вычислений
  в каждом процессе.

## 3. Описание алгоритма (последовательный)

Последовательная реализация основана на алгоритме LSD (Least Significant Digit):

1. Определение диапазона значений для учета отрицательных чисел.
2. Поразрядный проход (обычно по 10-ичному основанию или по байтам) с использованием сортировки подсчетом.
3. Обеспечение стабильности сортировки на каждом этапе.

## 4. Схема распараллеливания

Реализована гибридная схема "MPI + OpenMP":

1. **Распределение (MPI)**: Корневой процесс (rank 0) делит входной массив на части и рассылает их остальным
   процессам с помощью `MPI_Scatterv`.
2. **Локальная параллельная сортировка (OpenMP)**: Внутри каждого MPI-процесса полученный блок данных сортируется
   методом `LocalRadixSort`, который распараллелен с помощью OpenMP-директив (обработка разрядов и слияние
   подблоков).
3. **Сбор данных (MPI)**: Отсортированные блоки собираются на корневом процессе (`MPI_Gatherv`).
4. **Финальное слияние**: Корневой процесс выполняет итоговое слияние (merge) отсортированных частей в единый массив.

## 5. Детали реализации

- **Класс**: `VotincevDRadixMergeSortALL`.
- **MPI-функции**: `MPI_Bcast` для передачи размера, `MPI_Scatterv`/`MPI_Gatherv` для работы с блоками разного
  размера (если данные не делятся нацело).
- **OpenMP-ускорение**: Внутри `LocalRadixSort` используются `#pragma omp parallel` для разделения локального
  блока между потоками.
- **Память**: Каждый процесс выделяет память только под свою часть данных и временные буферы, что оптимизирует
  использование ОЗУ на распределенных системах.

## 6. Окружение

- **Оборудование**:
  - ЦП: Intel(R) Core(TM) i5-8400 CPU @ 2.80GHz (6 ядер)
  - ОЗУ: 16 ГБ
- **ПО**: Windows 10, MSVC, Microsoft MPI, OpenMP 2.0.
- **Тип сборки**: Release.
- **Данные**: Случайный вектор из 1 000 000 элементов `int32_t`.

## 7. Результаты

### 7.0 Воспроизводимость

Для запуска используются команды:

1) Конфигурация и билд
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
2) Func тесты
python scripts/run_tests.py --running-type=processes --counts 1 2 3 4
3) Perf тесты
$env:PPC_NUM_THREADS = "4"
$env:PPC_NUM_PROC = "2"
python scripts/run_tests.py --running-type="performance"

### 7.1 Корректность

Тестирование проводилось через Google Test в среде MPI. Проверялось выполнение на разном количестве процессов
и потоков. Результаты во всех случаях совпали с эталонной сортировкой `std::sort`.

### 7.2 Производительность

В таблице ниже представлены замеры времени выполнения для 1 000 000 элементов `int32_t`. Эксперименты проводились
с разным соотношением MPI-рангов и OpenMP-потоков для выявления наиболее эффективной конфигурации.

Количество элементов в массиве: 1 000 000

Время замеров - среднее. Режим task (не pipeline).

(эффективность считалась относительно общего числа ресурсов)

| Режим | Ranks (MPI) | Threads per Rank (OMP) | Всего потоков | Время, сек | Ускорение | Эффективность |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **seq** | 1 | 1 | 1 | 2.414 | 1.00 | 100.0% |
| **all** | 2 | 1 | 2 | 1.631 | **1.48** | 74.0% |
| **all** | 1 | 2 | 2 | 1.588 | **1.52** | 76.0% |
| **all** | 2 | 2 | 4 | 1.305 | 1.85 | 46.2% |
| **all** | 2 | 4 | 8 | 1.252 | 1.93 | 24.1% |

## 8. Заключение

Гибридная реализация MPI+OpenMP продемонстрировала эффективность, близкую к чистому OpenMP на одном узле.
Основное преимущество технологии ALL проявляется при масштабировании на несколько физических узлов кластера,
где MPI обеспечивает межсетевое взаимодействие, а OpenMP — эффективную работу с общей памятью внутри узла.

## 9. Список литературы

1. Сысоев А. В. Курс лекций "Параллельное программирование".
2. Официальная документация MPI

## Приложение

```cpp

// Если массив пустой - сортировать нечего
bool VotincevDRadixMergeSortALL::ValidationImpl() {
  return !GetInput().empty();
}

// Препроцессинга нет для сортировки
// (для изображений могло бы быть что-то вроде перевода в ч/б)
bool VotincevDRadixMergeSortALL::PreProcessingImpl() {
  return true;
}

// Постобработка не проводится
bool VotincevDRadixMergeSortALL::PostProcessingImpl() {
  return true;
}

bool VotincevDRadixMergeSortALL::RunImpl() {
  int32_t rank = 0;
  int32_t size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int32_t n = (rank == 0) ? static_cast<int32_t>(input_.size()) : 0;
  MPI_Bcast(&n, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);

  std::vector<int32_t> send_counts(static_cast<size_t>(size));
  std::vector<int32_t> displacements(static_cast<size_t>(size));
  int32_t items = n / size;
  int32_t rem = n % size;

  for (int32_t i = 0; i < size; ++i) {
    send_counts.at(static_cast<size_t>(i)) = items + (i < rem ? 1 : 0);
    displacements.at(static_cast<size_t>(i)) =
        (i == 0) ? 0 : displacements.at(static_cast<size_t>(i - 1)) + send_counts.at(static_cast<size_t>(i - 1));
  }

  auto local_n = send_counts.at(static_cast<size_t>(rank));
  std::vector<uint32_t> local_data(static_cast<size_t>(local_n));

  int32_t min_val = ScatterData(rank, n, local_n, send_counts, displacements, local_data);

  OmpLocalSortAndMerge(local_data);

  std::vector<uint32_t> gathered_data;
  if (rank == 0) {
    gathered_data.resize(static_cast<size_t>(n));
  }

  MPI_Gatherv(local_data.data(), local_n, MPI_UINT32_T, gathered_data.data(), send_counts.data(), displacements.data(),
              MPI_UINT32_T, 0, MPI_COMM_WORLD);

  FinalMergeAndFormat(rank, size, n, min_val, gathered_data, displacements);

  if (rank == 0) {
    GetOutput() = std::move(output_);
  }

  return true;
}
