#include "ivanova_p_marking_components_on_binary_image/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <string>
#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "ivanova_p_marking_components_on_binary_image/data/image_generator.hpp"

namespace ivanova_p_marking_components_on_binary_image {

IvanovaPMarkingComponentsOnBinaryImageTBB::IvanovaPMarkingComponentsOnBinaryImageTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();

  test_image.width = 0;
  test_image.height = 0;
  test_image.data.clear();
}

bool IvanovaPMarkingComponentsOnBinaryImageTBB::ValidationImpl() {
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
      test_image = CreateTestImage(500, 500, test_case);
    }
  }
  return test_image.width > 0 && !test_image.data.empty();
}

bool IvanovaPMarkingComponentsOnBinaryImageTBB::PreProcessingImpl() {
  input_image_ = test_image;
  width_ = input_image_.width;
  height_ = input_image_.height;
  int total_pixels = width_ * height_;

  labels_.assign(total_pixels, 0);

  // Инициализируем DSU: размер +1, так как метки теперь начинаются с 1 (idx + 1)
  parent_.resize(total_pixels + 1);
  for (int i = 0; i <= total_pixels; ++i) {
    parent_[i] = i;
  }

  current_label_ = 0;
  return true;
}

int IvanovaPMarkingComponentsOnBinaryImageTBB::FindRoot(int label) {
  int root = label;
  while (parent_[root] != root) {
    root = parent_[root];
  }

  // Сжатие путей (Path compression) для ускорения последующих поисков
  int current = label;
  while (parent_[current] != root) {
    int next = parent_[current];
    parent_[current] = root;
    current = next;
  }
  return root;
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::UnionLabels(int label1, int label2) {
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

void IvanovaPMarkingComponentsOnBinaryImageTBB::ProcessPixel(int xx, int yy, int idx) {
  int left_label = (xx > 0) ? labels_[idx - 1] : 0;
  int top_label = (yy > 0) ? labels_[idx - width_] : 0;

  bool left_exists = (left_label != 0);
  bool top_exists = (top_label != 0);

  if (!left_exists && !top_exists) {
    current_label_++;
    labels_[idx] = current_label_;
    parent_[current_label_] = current_label_;
  } else {
    int label = left_exists ? left_label : top_label;
    labels_[idx] = label;

    if (left_exists && top_exists && left_label != top_label) {
      UnionLabels(std::min(left_label, top_label), std::max(left_label, top_label));
    }
  }
}

int IvanovaPMarkingComponentsOnBinaryImageTBB::FindLocalRoot(int label) {
  int root = label;
  while (parent_[root] != root) {
    root = parent_[root];
  }
  return root;
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::UnionLocalRoots(int root1, int root2) {
  if (root1 != root2) {
    if (root1 < root2) {
      parent_[root2] = root1;
    } else {
      parent_[root1] = root2;
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::ProcessStripePixel(int xx, int yy, int idx, int start_row) {
  if (input_image_.data[idx] == 0) {
    return;
  }

  int left_label = (xx > 0) ? labels_[idx - 1] : 0;
  int top_label = (yy > start_row) ? labels_[idx - width_] : 0;

  if (left_label == 0 && top_label == 0) {
    labels_[idx] = idx + 1;
  } else {
    int label = (left_label != 0) ? left_label : top_label;
    labels_[idx] = label;

    if (left_label != 0 && top_label != 0 && left_label != top_label) {
      int root1 = FindLocalRoot(left_label);
      int root2 = FindLocalRoot(top_label);
      UnionLocalRoots(root1, root2);
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::FirstPass() {
  int num_threads = std::max(1, tbb::this_task_arena::max_concurrency());
  int rows_per_thread = (height_ + num_threads - 1) / num_threads;

  // Фаза 1: Истинная параллельная обработка независимых полос
  tbb::parallel_for(0, num_threads, [&](int thread_id) {
    int start_row = thread_id * rows_per_thread;
    int end_row = std::min(start_row + rows_per_thread, height_);
    if (start_row >= height_) {
      return;
    }

    for (int yy = start_row; yy < end_row; ++yy) {
      for (int xx = 0; xx < width_; ++xx) {
        int idx = (yy * width_) + xx;
        ProcessStripePixel(xx, yy, idx, start_row);
      }
    }
  });

  // Фаза 2: Быстрое последовательное объединение на стыках полос
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

  current_label_ = 1;
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::SecondPass() {
  int total_pixels = width_ * height_;

  // 1. Параллельно сжимаем пути для всех пикселей
  tbb::parallel_for(0, total_pixels, [&](int i) {
    if (labels_[i] != 0) {
      int root = labels_[i];
      while (parent_[root] != root) {
        root = parent_[root];
      }
      labels_[i] = root;
    }
  });

  // 2. Быстрый проход для создания непрерывных ID
  // ИСПРАВЛЕНИЕ: Вектор должен вмещать метки вплоть до total_pixels
  std::vector<int> root_mapping(total_pixels + 1, 0);
  int next_label = 1;
  for (int i = 0; i < total_pixels; ++i) {
    if (labels_[i] != 0) {
      int r = labels_[i];
      if (root_mapping[r] == 0) {
        root_mapping[r] = next_label++;
      }
    }
  }
  current_label_ = next_label - 1;

  // 3. Параллельно присваиваем новые нормализованные метки
  tbb::parallel_for(0, total_pixels, [&](int i) {
    if (labels_[i] != 0) {
      labels_[i] = root_mapping[labels_[i]];
    }
  });
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::InitLabelsTbb(int /*total_pixels*/) {
  // Не используется
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::MergeHorizontalPairsTbb() {
  // Не используется
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::MergeVerticalPairsTbb() {
  // Не используется
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::FinalizeRootsTbb(int /*total_pixels*/) {
  // Не используется
}

void IvanovaPMarkingComponentsOnBinaryImageTBB::NormalizeLabelsTbb(int /*total_pixels*/) {
  // Не используется
}

bool IvanovaPMarkingComponentsOnBinaryImageTBB::RunImpl() {
  if (width_ <= 0 || height_ <= 0) {
    return false;
  }

  FirstPass();

  if (current_label_ > 0) {
    SecondPass();
  }

  return true;
}

bool IvanovaPMarkingComponentsOnBinaryImageTBB::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();
  output.reserve(3 + labels_.size());
  output.push_back(width_);
  output.push_back(height_);
  output.push_back(current_label_);
  for (int l : labels_) {
    output.push_back(l);
  }
  return true;
}

}  // namespace ivanova_p_marking_components_on_binary_image
