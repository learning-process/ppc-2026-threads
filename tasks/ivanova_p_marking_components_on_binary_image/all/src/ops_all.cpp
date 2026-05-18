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
          break;
      }
      test_image = LoadImageFromTxt(filename);
    } else {
      int size = ExtractImageSize(test_case);
      test_image = CreateTestImage(size, size, test_case);
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
          break;
      }
      test_image = LoadImageFromTxt(filename);
    } else {
      int size = ExtractImageSize(test_case);
      test_image = CreateTestImage(size, size, test_case);
    }
  }

  input_image_ = test_image;
  width_ = input_image_.width;
  height_ = input_image_.height;

  int total_pixels = width_ * height_;
  labels_.assign(total_pixels, 0);

  // Выделяем память под глобальный DSU с запасом на все потоки
  int num_threads = tbb::this_task_arena::max_concurrency();
  parent_.resize((static_cast<size_t>(num_threads) * static_cast<size_t>(total_pixels)) + 1);
  for (size_t i = 0; i < parent_.size(); ++i) {
    parent_[i] = static_cast<int>(i);
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

int IvanovaPMarkingComponentsOnBinaryImageALL::FindLocalRootAll(int label, const std::vector<int> &local_parent) {
  int root = label;
  while (local_parent[static_cast<size_t>(root)] != root) {
    root = local_parent[static_cast<size_t>(root)];
  }
  return root;
}

void IvanovaPMarkingComponentsOnBinaryImageALL::UnionLocalLabelsAll(int label1, int label2,
                                                                    std::vector<int> &local_parent) {
  int root1 = FindLocalRootAll(label1, local_parent);
  int root2 = FindLocalRootAll(label2, local_parent);
  if (root1 != root2) {
    if (root1 < root2) {
      local_parent[static_cast<size_t>(root2)] = root1;
    } else {
      local_parent[static_cast<size_t>(root1)] = root2;
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::ProcessStripePixelAll(int xx, int yy, int idx, int start_row,
                                                                      std::vector<int> &local_parent,
                                                                      int &local_label) {
  if (input_image_.data[static_cast<size_t>(idx)] == 0) {
    return;
  }

  int left_label = (xx > 0) ? labels_[static_cast<size_t>(idx - 1)] : 0;
  int top_label = (yy > start_row) ? labels_[static_cast<size_t>(idx - width_)] : 0;

  bool left_exists = (left_label != 0);
  bool top_exists = (top_label != 0);

  if (!left_exists && !top_exists) {
    local_label++;
    labels_[static_cast<size_t>(idx)] = local_label;
  } else {
    int label = left_exists ? left_label : top_label;
    labels_[static_cast<size_t>(idx)] = label;

    if (left_exists && top_exists && left_label != top_label) {
      UnionLocalLabelsAll(left_label, top_label, local_parent);
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::InitializeLocalParent(std::vector<int> &local_parent, int max_labels) {
  local_parent.resize(static_cast<size_t>(max_labels));
  for (int i = 0; i < max_labels; ++i) {
    local_parent[static_cast<size_t>(i)] = i;
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::ProcessStripe(int start_row, int end_row,
                                                              std::vector<int> &local_parent, int &local_label) {
  for (int yy = start_row; yy < end_row; ++yy) {
    for (int xx = 0; xx < width_; ++xx) {
      int idx = (yy * width_) + xx;
      ProcessStripePixelAll(xx, yy, idx, start_row, local_parent, local_label);
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::MergeBoundaries(int num_threads, int rows_per_thread) {
  for (int thread_id = 0; thread_id < num_threads - 1; ++thread_id) {
    int boundary_row = (thread_id + 1) * rows_per_thread;
    if (boundary_row >= height_) {
      continue;
    }

    for (int xx = 0; xx < width_; ++xx) {
      int top_idx = ((boundary_row - 1) * width_) + xx;
      int bottom_idx = (boundary_row * width_) + xx;

      int top_label = labels_[static_cast<size_t>(top_idx)];
      int bottom_label = labels_[static_cast<size_t>(bottom_idx)];

      if (top_label != 0 && bottom_label != 0 && top_label != bottom_label) {
        UnionLabels(top_label, bottom_label);
      }
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::MergeLocalParents(const std::vector<std::vector<int>> &local_parents,
                                                                  const std::vector<int> &local_labels, int num_threads,
                                                                  int total_pixels) {
  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    if (local_parents[static_cast<size_t>(thread_id)].empty()) {
      continue;
    }

    int start_lbl = (thread_id * total_pixels) + 1;
    int end_lbl = local_labels[static_cast<size_t>(thread_id)];

    for (int label = start_lbl; label <= end_lbl; ++label) {
      int p = local_parents[static_cast<size_t>(thread_id)][static_cast<size_t>(label)];
      if (p != label) {
        UnionLabels(label, p);
      }
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::FirstPass() {
  int num_threads = tbb::this_task_arena::max_concurrency();
  int rows_per_thread = (height_ + num_threads - 1) / num_threads;
  int total_pixels = width_ * height_;

  std::vector<std::vector<int>> local_parents(static_cast<size_t>(num_threads));
  std::vector<int> local_labels(static_cast<size_t>(num_threads), 0);

  // Фаза 1: Параллельная обработка полос с TBB
  tbb::parallel_for(0, num_threads, [&](int thread_id) {
    int start_row = thread_id * rows_per_thread;
    int end_row = std::min(start_row + rows_per_thread, height_);
    if (start_row >= height_) {
      return;
    }

    int max_possible_labels = (num_threads * total_pixels) + 1;
    InitializeLocalParent(local_parents[static_cast<size_t>(thread_id)], max_possible_labels);

    int &local_label = local_labels[static_cast<size_t>(thread_id)];
    local_label = thread_id * total_pixels;

    ProcessStripe(start_row, end_row, local_parents[static_cast<size_t>(thread_id)], local_label);
  });

  // Фаза 2: Последовательное объединение границ между полосами
  MergeBoundaries(num_threads, rows_per_thread);

  // Фаза 3: Перенос локальных связей потоков в глобальный вектор DSU
  MergeLocalParents(local_parents, local_labels, num_threads, total_pixels);

  current_label_ = 0;
  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    current_label_ = std::max(current_label_, local_labels[static_cast<size_t>(thread_id)]);
  }
}

void IvanovaPMarkingComponentsOnBinaryImageALL::SecondPass() {
  int num_threads = tbb::this_task_arena::max_concurrency();
  int total_pixels = width_ * height_;

  // Вектор вместо unordered_map для нормализации меток
  std::vector<int> new_labels((static_cast<size_t>(num_threads) * static_cast<size_t>(total_pixels)) + 1, 0);
  int next_label = 1;

  for (int &label : labels_) {
    if (label != 0) {
      int root = FindRoot(label);

      if (new_labels[static_cast<size_t>(root)] == 0) {
        new_labels[static_cast<size_t>(root)] = next_label++;
      }

      label = new_labels[static_cast<size_t>(root)];
    }
  }

  current_label_ = next_label - 1;
}

void IvanovaPMarkingComponentsOnBinaryImageALL::InitLabelsAll(int /*unused*/) {}
void IvanovaPMarkingComponentsOnBinaryImageALL::MergeHorizontalPairsAll() {}
void IvanovaPMarkingComponentsOnBinaryImageALL::MergeVerticalPairsAll() {}
void IvanovaPMarkingComponentsOnBinaryImageALL::FinalizeRootsAll(int /*unused*/) {}
void IvanovaPMarkingComponentsOnBinaryImageALL::NormalizeLabelsAll(int /*unused*/) {}

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
