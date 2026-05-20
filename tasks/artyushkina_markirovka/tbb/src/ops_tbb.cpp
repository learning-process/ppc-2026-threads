#include "artyushkina_markirovka/tbb/include/ops_tbb.hpp"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include "artyushkina_markirovka/common/include/common.hpp"

namespace artyushkina_markirovka {

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

  int total_pixels = rows_ * cols_;
  labels_.assign(total_pixels, 0);
  parent_.resize(total_pixels + 1);
  for (int i = 0; i <= total_pixels; ++i) {
    parent_[i] = i;
  }

  current_label_ = 0;
  return true;
}

int MarkingComponentsTBB::FindRoot(int label) {
  // Итеративный поиск корня со сжатием пути
  int root = label;
  while (parent_[root] != root) {
    root = parent_[root];
  }

  // Сжатие пути
  int current = label;
  while (parent_[current] != current) {
    int next = parent_[current];
    parent_[current] = root;
    current = next;
  }
  return root;
}

void MarkingComponentsTBB::UnionLabels(int label1, int label2) {
  tbb::spin_mutex::scoped_lock lock(dsu_mutex_);
  int root1 = FindRoot(label1);
  int root2 = FindRoot(label2);
  if (root1 != root2) {
    // Объединяем по рангу (меньший корень становится родителем)
    if (root1 < root2) {
      parent_[root2] = root1;
    } else {
      parent_[root1] = root2;
    }
  }
}

void MarkingComponentsTBB::InitLabelsTbb() {
  int total_pixels = rows_ * cols_;
  tbb::parallel_for(0, total_pixels, [this](int idx) {
    size_t input_idx = static_cast<size_t>(idx) + 2;
    // 0 = объект, не-0 = фон
    if (input_[input_idx] == 0) {
      labels_[idx] = idx + 1;  // Временная метка
    }
  });
}

void MarkingComponentsTBB::MergeHorizontalPairsTbb() {
  // Для каждой строки обрабатываем горизонтальные связи
  tbb::parallel_for(0, rows_, [this](int y_coord) {
    for (int x_coord = 0; x_coord < cols_ - 1; ++x_coord) {
      int idx = (y_coord * cols_) + x_coord;
      if (labels_[idx] != 0 && labels_[idx + 1] != 0) {
        UnionLabels(labels_[idx], labels_[idx + 1]);
      }
    }
  });
}

void MarkingComponentsTBB::MergeVerticalPairsTbb() {
  // Для каждого столбца обрабатываем вертикальные связи
  tbb::parallel_for(0, cols_, [this](int x_coord) {
    for (int y_coord = 0; y_coord < rows_ - 1; ++y_coord) {
      int idx = (y_coord * cols_) + x_coord;
      if (labels_[idx] != 0 && labels_[idx + cols_] != 0) {
        UnionLabels(labels_[idx], labels_[idx + cols_]);
      }
    }
  });
}

void MarkingComponentsTBB::MergeDiagonalPairsTbb() {
  // Обрабатываем диагональные связи (8-связность)
  // Для верхней и левой диагоналей
  tbb::parallel_for(0, rows_, [this](int y_coord) {
    for (int x_coord = 0; x_coord < cols_; ++x_coord) {
      int idx = (y_coord * cols_) + x_coord;
      if (labels_[idx] == 0) {
        continue;
      }

      // Верхняя диагональ (северо-запад)
      if (y_coord > 0 && x_coord > 0 && labels_[idx - cols_ - 1] != 0) {
        UnionLabels(labels_[idx], labels_[idx - cols_ - 1]);
      }
      // Верхняя диагональ (северо-восток)
      if (y_coord > 0 && x_coord < cols_ - 1 && labels_[idx - cols_ + 1] != 0) {
        UnionLabels(labels_[idx], labels_[idx - cols_ + 1]);
      }
    }
  });
}

void MarkingComponentsTBB::FinalizeRootsTbb() {
  int total_pixels = rows_ * cols_;
  tbb::parallel_for(0, total_pixels, [this](int i) {
    if (labels_[i] != 0) {
      labels_[i] = FindRoot(labels_[i]);
    }
  });
}

void MarkingComponentsTBB::NormalizeLabelsTbb() {
  int total_pixels = rows_ * cols_;

  // Сначала собираем все уникальные корни
  std::vector<int> unique_roots;
  for (int i = 0; i < total_pixels; ++i) {
    if (labels_[i] != 0) {
      unique_roots.push_back(labels_[i]);
    }
  }

  // Сортируем и удаляем дубликаты
  std::sort(unique_roots.begin(), unique_roots.end());
  unique_roots.erase(std::unique(unique_roots.begin(), unique_roots.end()), unique_roots.end());

  // Создаем отображение
  std::vector<int> mapping(total_pixels + 1, 0);
  for (size_t i = 0; i < unique_roots.size(); ++i) {
    mapping[unique_roots[i]] = static_cast<int>(i + 1);
  }

  // Применяем отображение
  tbb::parallel_for(0, total_pixels, [this, &mapping](int i) {
    if (labels_[i] != 0) {
      labels_[i] = mapping[labels_[i]];
    }
  });

  current_label_ = static_cast<int>(unique_roots.size());
}

bool MarkingComponentsTBB::RunImpl() {
  int total_pixels = rows_ * cols_;
  if (total_pixels <= 0) {
    return true;
  }

  InitLabelsTbb();
  MergeHorizontalPairsTbb();
  MergeVerticalPairsTbb();
  MergeDiagonalPairsTbb();  // Добавляем обработку диагональных связей для 8-связности
  FinalizeRootsTbb();
  NormalizeLabelsTbb();

  return true;
}

bool MarkingComponentsTBB::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();

  output.push_back(static_cast<uint8_t>(rows_));
  output.push_back(static_cast<uint8_t>(cols_));

  for (int label : labels_) {
    output.push_back(static_cast<uint8_t>(label));
  }

  return true;
}

}  // namespace artyushkina_markirovka
