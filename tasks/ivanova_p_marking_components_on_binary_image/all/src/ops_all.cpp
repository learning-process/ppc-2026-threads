#include "ivanova_p_marking_components_on_binary_image/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "ivanova_p_marking_components_on_binary_image/data/image_generator.hpp"

namespace ivanova_p_marking_components_on_binary_image {

IvanovaPMarkingComponentsOnBinaryImageALL::IvanovaPMarkingComponentsOnBinaryImageALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();

  test_image.width = 0;
  test_image.height = 0;
  test_image.data.clear();
}

bool IvanovaPMarkingComponentsOnBinaryImageALL::ValidationImpl() {
  if (test_image.width <= 0 || test_image.height <= 0) {
    int test_case = GetInput();

    if (test_case >= 11 && test_case <= 14) {
      std::string filename;
      switch (test_case) {
        case 11:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image.txt";
          break;
        case 12:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image2.txt";
          break;
        case 13:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image3.txt";
          break;
        case 14:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image4.txt";
          break;
        default:
          filename = "";
      }
      test_image = LoadImageFromTxt(filename);
    } else {
      const int width = 500;
      const int height = 500;
      test_image = CreateTestImage(width, height, test_case);
    }
  }

  if (test_image.width <= 0 || test_image.height <= 0 || test_image.data.empty()) {
    return false;
  }
  if (test_image.data.size() != static_cast<size_t>(test_image.width) * static_cast<size_t>(test_image.height)) {
    return false;
  }
  return true;
}

bool IvanovaPMarkingComponentsOnBinaryImageALL::PreProcessingImpl() {
  if (test_image.width <= 0 || test_image.height <= 0) {
    // Дублирующий блок из оригинального кода оставлен для сохранения логики
    int test_case = GetInput();
    if (test_case >= 11 && test_case <= 14) {
      std::string filename;
      switch (test_case) {
        case 11:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image.txt";
          break;
        case 12:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image2.txt";
          break;
        case 13:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image3.txt";
          break;
        case 14:
          filename = "tasks/ivanova_p_marking_components_on_binary_image/data/image4.txt";
          break;
        default:
          filename = "";
      }
      test_image = LoadImageFromTxt(filename);
    } else {
      const int width = 500;
      const int height = 500;
      test_image = CreateTestImage(width, height, test_case);
    }
  }

  input_image_ = test_image;
  width_ = input_image_.width;
  height_ = input_image_.height;

  int total_pixels = width_ * height_;
  labels_.assign(total_pixels, 0);

  // Выделяем память под глобальный DSU с запасом на все потоки
  int num_threads = tbb::this_task_arena::max_concurrency();
  parent_.resize(num_threads * total_pixels + 1);
  for (size_t i = 0; i < parent_.size(); ++i) {
    parent_[i] = i;
  }

  current_label_ = 0;
  return true;
}

int IvanovaPMarkingComponentsOnBinaryImageALL::FindRoot(int label) {
  int root = label;
  while (parent_[root] != root) {
    root = parent_[root];
  }

  // Сжатие путей (Path compression)
  int current = label;
  while (parent_[current] != root) {
    int next = parent_[current];
    parent_[current] = root;
    current = next;
  }
  return root;
}

void IvanovaPMarkingComponentsOnBinaryImageALL::UnionLabels(int label1, int label2) {
  int root1 = FindRoot(label1);
  int root2 = FindRoot(label2);

  if (root1 != root2) {
    if (root1 < root2) {
      parent_[root2] = root1;
    } else {
      parent_[root1] = root2;
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::ProcessPixel(int /*xx*/, int /*yy*/, int /*idx*/) {
  // Не используется в основном параллельном проходе, оставлено для интерфейса класса
}

void IvanovaPMarkingComponentsOnBinaryImageALL::FirstPass() {
  int num_threads = tbb::this_task_arena::max_concurrency();
  int rows_per_thread = (height_ + num_threads - 1) / num_threads;
  int total_pixels = width_ * height_;

  // Локальные DSU для каждого потока теперь на векторах
  std::vector<std::vector<int>> local_parents(num_threads);
  std::vector<int> local_labels(num_threads, 0);

  // Фаза 1: Параллельная обработка полос с TBB
  tbb::parallel_for(0, num_threads, [&](int thread_id) {
    int start_row = thread_id * rows_per_thread;
    int end_row = std::min(start_row + rows_per_thread, height_);
    if (start_row >= height_) {
      return;
    }

    // Инициализация локального DSU вектора для конкретного потока
    int max_possible_labels = num_threads * total_pixels + 1;
    local_parents[thread_id].resize(max_possible_labels);
    for (int i = 0; i < max_possible_labels; ++i) {
      local_parents[thread_id][i] = i;
    }

    int &local_label = local_labels[thread_id];
    local_label = thread_id * total_pixels;  // Гарантия уникального диапазона меток

    for (int yy = start_row; yy < end_row; ++yy) {
      for (int xx = 0; xx < width_; ++xx) {
        int idx = (yy * width_) + xx;

        if (input_image_.data[idx] == 0) {
          continue;
        }

        int left_label = (xx > 0) ? labels_[idx - 1] : 0;
        int top_label = (yy > start_row) ? labels_[idx - width_] : 0;

        bool left_exists = (left_label != 0);
        bool top_exists = (top_label != 0);

        if (!left_exists && !top_exists) {
          local_label++;
          labels_[idx] = local_label;
        } else {
          int label = left_exists ? left_label : top_label;
          labels_[idx] = label;

          if (left_exists && top_exists && left_label != top_label) {
            auto &lp = local_parents[thread_id];
            int root1 = left_label;
            int root2 = top_label;
            while (lp[root1] != root1) {
              root1 = lp[root1];
            }
            while (lp[root2] != root2) {
              root2 = lp[root2];
            }
            if (root1 != root2) {
              if (root1 < root2) {
                lp[root2] = root1;
              } else {
                lp[root1] = root2;
              }
            }
          }
        }
      }
    }
  });

  // Фаза 2: Последовательное объединение границ между полосами на уровне процесса
  for (int thread_id = 0; thread_id < num_threads - 1; ++thread_id) {
    int boundary_row = (thread_id + 1) * rows_per_thread;
    if (boundary_row >= height_) {
      continue;
    }

    for (int xx = 0; xx < width_; ++xx) {
      int top_idx = (boundary_row - 1) * width_ + xx;
      int bottom_idx = boundary_row * width_ + xx;

      int top_label = labels_[top_idx];
      int bottom_label = labels_[bottom_idx];

      if (top_label != 0 && bottom_label != 0 && top_label != bottom_label) {
        UnionLabels(top_label, bottom_label);
      }
    }
  }

  // Фаза 3: Перенос локальных связей потоков в глобальный вектор DSU
  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    if (local_parents[thread_id].empty()) {
      continue;
    }

    int start_lbl = thread_id * total_pixels + 1;
    int end_lbl = local_labels[thread_id];

    for (int label = start_lbl; label <= end_lbl; ++label) {
      int p = local_parents[thread_id][label];
      if (p != label) {
        UnionLabels(label, p);
      }
    }
  }

  current_label_ = 0;
  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    current_label_ = std::max(current_label_, local_labels[thread_id]);
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::SecondPass() {
  int num_threads = tbb::this_task_arena::max_concurrency();
  int total_pixels = width_ * height_;

  // Вектор вместо unordered_map для нормализации меток
  std::vector<int> new_labels(num_threads * total_pixels + 1, 0);
  int next_label = 1;

  for (int &label : labels_) {
    if (label != 0) {
      int root = FindRoot(label);

      if (new_labels[root] == 0) {
        new_labels[root] = next_label++;
      }

      label = new_labels[root];
    }
  }

  current_label_ = next_label - 1;
}

void IvanovaPMarkingComponentsOnBinaryImageALL::InitLabelsAll(int) {}
void IvanovaPMarkingComponentsOnBinaryImageALL::MergeHorizontalPairsAll() {}
void IvanovaPMarkingComponentsOnBinaryImageALL::MergeVerticalPairsAll() {}
void IvanovaPMarkingComponentsOnBinaryImageALL::FinalizeRootsAll(int) {}
void IvanovaPMarkingComponentsOnBinaryImageALL::NormalizeLabelsAll(int) {}

bool IvanovaPMarkingComponentsOnBinaryImageALL::RunImpl() {
  int total_pixels = width_ * height_;
  if (total_pixels <= 0) {
    return false;
  }

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  FirstPass();

  if (current_label_ > 0) {
    SecondPass();
  }

  if (world_size > 1) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&current_label_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(labels_.data(), static_cast<int>(labels_.size()), MPI_INT, 0, MPI_COMM_WORLD);
  }

  return true;
}

bool IvanovaPMarkingComponentsOnBinaryImageALL::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();

  output.push_back(width_);
  output.push_back(height_);
  output.push_back(current_label_);

  for (int label : labels_) {
    output.push_back(label);
  }

  return true;
}

}  // namespace ivanova_p_marking_components_on_binary_image
