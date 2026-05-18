#include "ivanova_p_marking_components_on_binary_image/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <thread>
#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "ivanova_p_marking_components_on_binary_image/data/image_generator.hpp"
#include "util/include/util.hpp"

namespace ivanova_p_marking_components_on_binary_image {

IvanovaPMarkingComponentsOnBinaryImageSTL::IvanovaPMarkingComponentsOnBinaryImageSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();

  test_image.width = 0;
  test_image.height = 0;
  test_image.data.clear();
}

bool IvanovaPMarkingComponentsOnBinaryImageSTL::ValidationImpl() {
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

bool IvanovaPMarkingComponentsOnBinaryImageSTL::PreProcessingImpl() {
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

  // Инициализация глобального вектора DSU
  int num_threads = ppc::util::GetNumThreads();
  parent_.resize((static_cast<size_t>(num_threads) * static_cast<size_t>(total_pixels)) + 1);
  for (size_t i = 0; i < parent_.size(); ++i) {
    parent_[i] = static_cast<int>(i);
  }

  current_label_ = 0;
  return true;
}

int IvanovaPMarkingComponentsOnBinaryImageSTL::FindRoot(int label) {
  int root = label;
  while (parent_[root] != root) {
    root = parent_[root];
  }

  int current = label;
  while (parent_[current] != root) {
    int next = parent_[current];
    parent_[current] = root;
    current = next;
  }
  return root;
}

void IvanovaPMarkingComponentsOnBinaryImageSTL::UnionLabels(int label1, int label2) {
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

void IvanovaPMarkingComponentsOnBinaryImageSTL::ProcessPixel(int /*xx*/, int /*yy*/, int /*idx*/) {}

int IvanovaPMarkingComponentsOnBinaryImageSTL::FindLocalRootStl(int label, const std::vector<int> &local_parent) {
  int root = label;
  while (local_parent[root] != root) {
    root = local_parent[root];
  }
  return root;
}

void IvanovaPMarkingComponentsOnBinaryImageSTL::UnionLocalLabelsStl(int label1, int label2,
                                                                    std::vector<int> &local_parent) {
  int root1 = FindLocalRootStl(label1, local_parent);
  int root2 = FindLocalRootStl(label2, local_parent);
  if (root1 != root2) {
    if (root1 < root2) {
      local_parent[root2] = root1;
    } else {
      local_parent[root1] = root2;
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageSTL::ProcessStripePixelStl(int xx, int yy, int idx, int start_row,
                                                                      std::vector<int> &local_parent,
                                                                      int &local_label) {
  if (input_image_.data[idx] == 0) {
    return;
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
      UnionLocalLabelsStl(left_label, top_label, local_parent);
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageSTL::FirstPass() {
  int num_threads = ppc::util::GetNumThreads();
  int rows_per_thread = (height_ + num_threads - 1) / num_threads;
  int total_pixels = width_ * height_;

  // Локальные DSU переведены на векторы
  std::vector<std::vector<int>> local_parents(num_threads);
  std::vector<int> local_labels(num_threads, 0);
  std::vector<std::thread> threads;

  threads.reserve(static_cast<size_t>(num_threads));
  // Фаза 1: Параллельная обработка полос через стандартные потоки
  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    threads.emplace_back(
        [this, thread_id, rows_per_thread, &local_parents, &local_labels, num_threads, total_pixels]() {
      int start_row = thread_id * rows_per_thread;
      int end_row = std::min(start_row + rows_per_thread, height_);
      if (start_row >= height_) {
        return;
      }

      int max_possible_labels = (num_threads * total_pixels) + 1;
      local_parents[thread_id].resize(static_cast<size_t>(max_possible_labels));
      for (int i = 0; i < max_possible_labels; ++i) {
        local_parents[thread_id][static_cast<size_t>(i)] = i;
      }

      int &local_label = local_labels[static_cast<size_t>(thread_id)];
      local_label = thread_id * total_pixels;

      for (int yy = start_row; yy < end_row; ++yy) {
        for (int xx = 0; xx < width_; ++xx) {
          int idx = (yy * width_) + xx;
          ProcessStripePixelStl(xx, yy, idx, start_row, local_parents[thread_id], local_label);
        }
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  // Фаза 2: Последовательное объединение смежных границ между полосами
  for (int thread_id = 0; thread_id < num_threads - 1; ++thread_id) {
    int boundary_row = (thread_id + 1) * rows_per_thread;
    if (boundary_row >= height_) {
      continue;
    }

    for (int xx = 0; xx < width_; ++xx) {
      int top_idx = ((boundary_row - 1) * width_) + xx;
      int bottom_idx = (boundary_row * width_) + xx;

      int top_label = labels_[top_idx];
      int bottom_label = labels_[bottom_idx];

      if (top_label != 0 && bottom_label != 0 && top_label != bottom_label) {
        UnionLabels(top_label, bottom_label);
      }
    }
  }

  // Фаза 3: Перенос локальных отношений эквивалентности в глобальный вектор
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

  current_label_ = 0;
  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    current_label_ = std::max(current_label_, local_labels[static_cast<size_t>(thread_id)]);
  }
}

void IvanovaPMarkingComponentsOnBinaryImageSTL::SecondPass() {
  int num_threads = ppc::util::GetNumThreads();
  int total_pixels = width_ * height_;

  // Убран локальный unordered_map, заменен плоским вектором разметки
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

void IvanovaPMarkingComponentsOnBinaryImageSTL::InitLabelsStl(int /*unused*/, int /*unused*/) {}
void IvanovaPMarkingComponentsOnBinaryImageSTL::MergeHorizontalPairsStl(int /*unused*/) {}
void IvanovaPMarkingComponentsOnBinaryImageSTL::MergeVerticalPairsStl(int /*unused*/) {}
void IvanovaPMarkingComponentsOnBinaryImageSTL::FinalizeRootsStl(int /*unused*/, int /*unused*/) {}
void IvanovaPMarkingComponentsOnBinaryImageSTL::NormalizeLabelsStl(int /*unused*/) {}

bool IvanovaPMarkingComponentsOnBinaryImageSTL::RunImpl() {
  int total_pixels = width_ * height_;
  if (total_pixels <= 0) {
    return false;
  }

  FirstPass();

  if (current_label_ > 0) {
    SecondPass();
  }

  return true;
}

bool IvanovaPMarkingComponentsOnBinaryImageSTL::PostProcessingImpl() {
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
