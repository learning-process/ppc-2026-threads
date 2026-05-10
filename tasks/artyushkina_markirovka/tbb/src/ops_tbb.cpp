#include "artyushkina_markirovka/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/mutex.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <map>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"

namespace artyushkina_markirovka {
namespace {

tbb::mutex union_mutex;

// Вспомогательная функция для проверки и добавления соседа
void AddNeighborIfValid(int neighbor_label, std::vector<int> &neighbor_labels) {
  if (neighbor_label != 0) {
    neighbor_labels.push_back(neighbor_label);
  }
}

// 8-связность: сбор меток соседей (только для объектов)
void CollectNeighborsLabels(int i, int j, const std::vector<std::vector<int>> &temp_labels,
                            std::vector<int> &neighbor_labels, int /*rows*/, int cols) {
  // Сосед сверху-слева (диагональ)
  if (i > 0 && j > 0) {
    AddNeighborIfValid(temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j - 1)], neighbor_labels);
  }
  // Сосед сверху
  if (i > 0) {
    AddNeighborIfValid(temp_labels[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)], neighbor_labels);
  }
  // Сосед сверху-справа (диагональ)
  if (i > 0 && j + 1 < cols) {
    auto row_idx = static_cast<std::size_t>(i - 1);
    auto col_idx = static_cast<std::size_t>(j + 1);
    AddNeighborIfValid(temp_labels[row_idx][col_idx], neighbor_labels);
  }
  // Сосед слева
  if (j > 0) {
    AddNeighborIfValid(temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)], neighbor_labels);
  }
}

int FindMinLabel(const std::vector<int> &labels) {
  if (labels.empty()) {
    return 0;
  }
  int min_label = labels[0];
  for (std::size_t k = 1; k < labels.size(); ++k) {
    min_label = std::min(min_label, labels[k]);
  }
  return min_label;
}

// Вынесем логику обработки пикселя в отдельную функцию для снижения когнитивной сложности
void ProcessPixel(int i, int j, const InType &input, int cols, std::vector<std::vector<int>> &temp_labels,
                  std::vector<int> &parent, std::atomic<int> &next_label) {
  auto idx = (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(j) + 2;

  // Проверка: 0 означает объект (черный), ненулевое - фон (белый)
  if (input[idx] != 0) {
    return;  // пропускаем фон
  }

  std::vector<int> neighbor_labels;
  neighbor_labels.reserve(4);

  // Собираем метки соседей с 8-связностью
  CollectNeighborsLabels(i, j, temp_labels, neighbor_labels, 0, cols);

  if (neighbor_labels.empty()) {
    // Нет соседей - новый компонент
    int label = next_label++;
    temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = label;

    tbb::mutex::scoped_lock lock(union_mutex);
    if (static_cast<std::size_t>(label) >= parent.size()) {
      parent.resize(static_cast<std::size_t>(label) + 1);
    }
    parent[static_cast<std::size_t>(label)] = label;
  } else {
    // Есть соседи - берем минимальную метку
    int min_label = FindMinLabel(neighbor_labels);
    temp_labels[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = min_label;

    // Объединяем все метки соседей
    for (int label : neighbor_labels) {
      if (label != min_label) {
        MarkingComponentsTBB::UnionLabels(parent, min_label, label);
      }
    }
  }
}

}  // namespace

MarkingComponentsTBB::MarkingComponentsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MarkingComponentsTBB::ValidationImpl() {
  return GetInput().size() >= 2;
}

bool MarkingComponentsTBB::PreProcessingImpl() {
  const auto &input = GetInput();
  rows_ = static_cast<int>(input[0]);
  cols_ = static_cast<int>(input[1]);
  input_ = input;

  labels_.clear();
  labels_.resize(static_cast<std::size_t>(rows_));
  for (int i = 0; i < rows_; ++i) {
    labels_[static_cast<std::size_t>(i)].assign(static_cast<std::size_t>(cols_), 0);
  }

  temp_labels_.clear();
  temp_labels_.resize(static_cast<std::size_t>(rows_));
  for (int i = 0; i < rows_; ++i) {
    temp_labels_[static_cast<std::size_t>(i)].assign(static_cast<std::size_t>(cols_), 0);
  }

  parent_.clear();
  parent_.push_back(0);
  next_label_ = 1;

  return true;
}

int MarkingComponentsTBB::FindRoot(std::vector<int> &parent, int label) {
  int current_label = label;
  while (parent[static_cast<std::size_t>(current_label)] != current_label) {
    parent[static_cast<std::size_t>(current_label)] =
        parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(current_label)])];
    current_label = parent[static_cast<std::size_t>(current_label)];
  }
  return current_label;
}

void MarkingComponentsTBB::UnionLabels(std::vector<int> &parent, int label1, int label2) {
  if (label1 == label2) {
    return;
  }

  tbb::mutex::scoped_lock lock(union_mutex);
  int root1 = FindRoot(parent, label1);
  int root2 = FindRoot(parent, label2);

  if (root1 != root2) {
    if (root1 < root2) {
      parent[static_cast<std::size_t>(root2)] = root1;
    } else {
      parent[static_cast<std::size_t>(root1)] = root2;
    }
  }
}

void MarkingComponentsTBB::ProcessFirstPass() {
  tbb::parallel_for(0, rows_, [&](int i) {
    for (int j = 0; j < cols_; ++j) {
      ProcessPixel(i, j, input_, cols_, temp_labels_, parent_, next_label_);
    }
  });
}

void MarkingComponentsTBB::ResolveEquivalences() {
  tbb::parallel_for(0, rows_, [&](int i) {
    for (int j = 0; j < cols_; ++j) {
      int &label = temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      if (label != 0) {
        label = FindRoot(parent_, label);
      }
    }
  });
}

void MarkingComponentsTBB::RemapLabels() {
  // Собираем уникальные метки
  std::vector<int> unique_labels;
  unique_labels.reserve(static_cast<std::size_t>(rows_) * static_cast<std::size_t>(cols_));

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      int label = temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      if (label != 0) {
        unique_labels.push_back(label);
      }
    }
  }

  std::sort(unique_labels.begin(), unique_labels.end());
  auto last = std::unique(unique_labels.begin(), unique_labels.end());
  unique_labels.erase(last, unique_labels.end());

  // Создаем отображение старых меток на новые (1, 2, 3, ...)
  std::map<int, int> label_mapping;
  int current_label = 1;
  for (int label : unique_labels) {
    label_mapping[label] = current_label++;
  }

  // Применяем отображение
  tbb::parallel_for(0, rows_, [&](int i) {
    for (int j = 0; j < cols_; ++j) {
      int label = temp_labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      if (label != 0) {
        labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = label_mapping[label];
      } else {
        labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = 0;
      }
    }
  });
}

bool MarkingComponentsTBB::RunImpl() {
  if (input_.size() < 2 || rows_ == 0 || cols_ == 0) {
    return false;
  }

  ProcessFirstPass();
  ResolveEquivalences();
  RemapLabels();

  return true;
}

bool MarkingComponentsTBB::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();

  output.push_back(rows_);
  output.push_back(cols_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      output.push_back(labels_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
    }
  }

  return true;
}

}  // namespace artyushkina_markirovka
