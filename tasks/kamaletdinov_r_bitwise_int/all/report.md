# Поразрядная сортировка для целых чисел — ALL

- Студент: Камалетдинов Рамзан Рамилевич
- Группа: 3823Б1ПР4
- Технология: ALL (MPI + OMP)
- Вариант: 17

## 1. Контекст

Гибридная версия сочетает межпроцессное взаимодействие через MPI и внутрипроцессное распараллеливание через OpenMP. Это
позволяет использовать ресурсы кластера или многопроцессорных систем.

## 2. Постановка задачи

Распределение данных между процессами MPI, локальная сортировка с использованием OpenMP и последующее слияние
результатов. Реализация основана на последовательной версии алгоритма.

## 3. Базовый алгоритм

Гибридный алгоритм реализует двухуровневую схему параллелизма, значительно усложняя структуру последовательного
алгоритма.

В отличие от SEQ-версии, где данные обрабатываются локально и последовательно:

- **Уровень 1 (Межпроцессный)**: Данные не просто лежат в памяти, а активно распределяются. Root-процесс разбивает
  массив на блоки и рассылает их через `MPI_Scatterv`. В последовательной версии этот этап отсутствует вовсе.
- **Уровень 2 (Внутрипроцессный)**: Каждый MPI-процесс на своем локальном блоке данных запускает поразрядную
  сортировку. Внутри локальной сортировки этапы подсчета гистограмм и перераспределения данных ускоряются с помощью
  OpenMP. Это превращает линейный проход SEQ в многопоточный.
- **Слияние**: После сбора всех отсортированных частей через `MPI_Gatherv`, root-процесс должен выполнить процедуру
  многопутевого слияния, так как локально отсортированные куски не образуют глобально отсортированный массив
  автоматически. В последовательной версии массив всегда остается единым целым.

## 4. Межпроцессная схема

- **Роли**: Rank 0 — координатор (рассылка, прием, слияние). Остальные — вычислители.
- **MPI-вызовы**: `MPI_Scatterv` для рассылки блоков разного размера (если $N$ не делится нацело), `MPI_Gatherv` для
  сбора отсортированных блоков.
- **Синхронизация**: Использование `MPI_Barrier` для замера чистого времени работы алгоритма.

## 5. Внутрипроцессная схема

Внутри каждого процесса используется OpenMP для распараллеливания поразрядной сортировки. Это позволяет эффективно
задействовать все ядра процессора на каждом узле.

## 6. Детали реализации

Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`
Основная сложность заключается в корректном расчете смещений для `MPI_Scatterv`/`MPI_Gatherv` и реализации
эффективного слияния на стороне root-процесса.

## 7. Проверка корректности

Сравнение финального массива на root-процессе с последовательным эталоном. Тестирование при различных конфигурациях
`PPC_NUM_PROC` и `PPC_NUM_THREADS`.

## 8. Экспериментальная среда

- **CPU**: Apple M1 (8 ядер)
- **Компилятор**: GCC 13.3.0
- **Build Type**: Release

## 9. Результаты

Замеры производительности:

**task_run**:

| Mode | Processes | Threads | Time, s  | Speedup | Efficiency |
|------|-----------|---------|----------|---------|------------|
| seq  | 1         | 1       | 0.087591 | 1.000   | N/A        |
| all  | 1         | 1       | 0.054988 | 1.593   | 159.3%     |
| all  | 1         | 3       | 0.064376 | 1.360   | 45.3%      |
| all  | 1         | 6       | 0.055734 | 1.572   | 26.2%      |
| all  | 1         | 9       | 0.055940 | 1.566   | 17.4%      |
| all  | 3         | 1       | 0.092298 | 0.949   | 31.6%      |
| all  | 3         | 3       | 0.092076 | 0.951   | 10.6%      |
| all  | 3         | 6       | 0.092941 | 0.942   | 5.2%       |
| all  | 3         | 9       | 0.093175 | 0.940   | 3.5%       |
| all  | 6         | 1       | 0.159347 | 0.549   | 9.2%       |
| all  | 6         | 3       | 0.162089 | 0.540   | 3.0%       |
| all  | 6         | 6       | 0.159670 | 0.548   | 1.5%       |
| all  | 6         | 9       | 0.160124 | 0.547   | 1.0%       |
| all  | 9         | 1       | 0.219271 | 0.399   | 4.4%       |
| all  | 9         | 3       | 0.223071 | 0.392   | 1.5%       |
| all  | 9         | 6       | 0.219047 | 0.400   | 0.7%       |
| all  | 9         | 9       | 0.228623 | 0.383   | 0.5%       |

**task_pipeline**:

| Mode | Processes | Threads | Time, s  | Speedup | Efficiency |
|------|-----------|---------|----------|---------|------------|
| seq  | 1         | 1       | 0.104688 | 1.000   | N/A        |
| all  | 1         | 1       | 0.057609 | 1.817   | 181.7%     |
| all  | 1         | 3       | 0.057171 | 1.831   | 61.0%      |
| all  | 1         | 6       | 0.056849 | 1.841   | 30.7%      |
| all  | 1         | 9       | 0.057336 | 1.825   | 20.3%      |
| all  | 3         | 1       | 0.094143 | 1.112   | 37.1%      |
| all  | 3         | 3       | 0.093321 | 1.121   | 12.5%      |
| all  | 3         | 6       | 0.095417 | 1.097   | 6.1%       |
| all  | 3         | 9       | 0.093692 | 1.117   | 4.1%       |
| all  | 6         | 1       | 0.155531 | 0.673   | 11.2%      |
| all  | 6         | 3       | 0.155829 | 0.671   | 3.7%       |
| all  | 6         | 6       | 0.155843 | 0.671   | 1.9%       |
| all  | 6         | 9       | 0.155828 | 0.671   | 1.2%       |
| all  | 9         | 1       | 0.224792 | 0.465   | 5.2%       |
| all  | 9         | 3       | 0.228265 | 0.458   | 1.7%       |
| all  | 9         | 6       | 0.225249 | 0.464   | 0.9%       |
| all  | 9         | 9       | 0.225761 | 0.463   | 0.6%       |

## 10. Выводы

Гибридная схема на данном объеме данных ($N$) не дает преимущества над чистыми потоковыми реализациями. Это
объясняется высокой стоимостью передачи данных между процессами и накладными расходами на финальное слияние. MPI
версия целесообразна только на очень больших объемах данных, превышающих объем оперативной памяти одного узла.

## 11. Код реализации

```cpp
void BitwiseSortALL(std::vector<int> &arr) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int total = static_cast<int>(arr.size());
  MPI_Bcast(&total, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total <= 1) return;

  std::vector<int> send_counts(size);
  std::vector<int> displs(size);
  for (int i = 0; i < size; i++) {
    send_counts[i] = (total / size) + (i < total % size ? 1 : 0);
    displs[i] = (i == 0) ? 0 : displs[i - 1] + send_counts[i - 1];
  }

  std::vector<int> local_data(send_counts[rank]);
  MPI_Scatterv(arr.data(), send_counts.data(), displs.data(), MPI_INT, 
               local_data.data(), send_counts[rank], MPI_INT, 0, MPI_COMM_WORLD);

  LocalBitwiseSort(local_data);

  std::vector<int> recv_counts(size);
  MPI_Gather(&send_counts[rank], 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::vector<int> all_data(total);
    std::vector<int> recv_displs(size);
    for (int i = 1; i < size; i++) {
      recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
    }
    MPI_Gatherv(local_data.data(), send_counts[rank], MPI_INT, all_data.data(), 
                recv_counts.data(), recv_displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> merged = {all_data.begin(), all_data.begin() + recv_counts[0]};
    for (int i = 1; i < size; i++) {
      std::vector<int> chunk(all_data.begin() + recv_displs[i], 
                             all_data.begin() + recv_displs[i] + recv_counts[i]);
      merged = MergeSorted(merged, chunk);
    }
    arr = merged;
  } else {
    MPI_Gatherv(local_data.data(), send_counts[rank], MPI_INT, nullptr, nullptr, 
                nullptr, MPI_INT, 0, MPI_COMM_WORLD);
    arr.resize(total);
  }

  MPI_Bcast(arr.data(), total, MPI_INT, 0, MPI_COMM_WORLD);
}

void CountingSortByDigitALL(std::vector<int> &arr, int exp) {
  const int n = static_cast<int>(arr.size());
  const int thread_count = omp_get_max_threads();
  std::vector<std::array<int, 10>> local_counts(thread_count);

#pragma omp parallel default(none) shared(arr, exp, local_counts, n) num_threads(thread_count)
  {
    const int tid = omp_get_thread_num();
    auto &current = local_counts.at(tid);
    current.fill(0);

#pragma omp for schedule(static)
    for (int i = 0; i < n; i++) {
      const int digit = (arr.at(i) / exp) % 10;
      current.at(digit)++;
    }
  }

  // сборка глобальной гистограммы и расчет смещений...
}
```
