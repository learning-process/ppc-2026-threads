#include "ivanova_p_marking_components_on_binary_image/seq/include/ops_seq.hpp"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

namespace ivanova_p_marking_components_on_binary_image {

IvanovaPMarkingComponentsOnBinaryImageSEQ::IvanovaPMarkingComponentsOnBinaryImageSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();

  // Очищаем глобальное изображение для нового теста
  test_image.width = 0;
  test_image.height = 0;
  test_image.data.clear();
}

bool IvanovaPMarkingComponentsOnBinaryImageSEQ::ValidationImpl() {
  // Если изображение еще не инициализировано, инициализируем его
  if (test_image.width <= 0 || test_image.height <= 0) {
    int test_case = GetInput();

    if (test_case >= 11 && test_case <= 14) {
      // Загружаем из файла
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
      // Создаем программно
      const int width = 100;
      const int height = 100;
      test_image = CreateTestImage(width, height, test_case);
    }
  }

  if (test_image.width <= 0 || test_image.height <= 0) {
    return false;
  }
  if (test_image.data.empty()) {
    return false;
  }
  if (test_image.data.size() != static_cast<size_t>(test_image.width * test_image.height)) {
    return false;
  }
  return true;
}

bool IvanovaPMarkingComponentsOnBinaryImageSEQ::PreProcessingImpl() {
  // Если изображение еще не инициализировано, инициализируем его
  if (test_image.width <= 0 || test_image.height <= 0) {
    // Используем входные данные для определения типа теста
    int test_case = GetInput();

    if (test_case >= 11 && test_case <= 14) {
      // Загружаем из файла
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
      // Создаем программно
      const int width = 100;
      const int height = 100;
      test_image = CreateTestImage(width, height, test_case);
    }
  }

  input_image_ = test_image;
  width_ = input_image_.width;
  height_ = input_image_.height;

  labels_.assign(static_cast<size_t>(width_ * height_), 0);
  parent_.clear();
  current_label_ = 0;

  return true;
}

int IvanovaPMarkingComponentsOnBinaryImageSEQ::FindRoot(int label) {
  if (!parent_.contains(label)) {
    parent_[label] = label;
    return label;
  }

  if (parent_[label] != label) {
    parent_[label] = FindRoot(parent_[label]);
  }

  return parent_[label];
}

void IvanovaPMarkingComponentsOnBinaryImageSEQ::UnionLabels(int label1, int label2) {
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

void IvanovaPMarkingComponentsOnBinaryImageSEQ::FirstPass() {
  for (int yy = 0; yy < height_; ++yy) {
    for (int xx = 0; xx < width_; ++xx) {
      int idx = yy * width_ + xx;

      if (input_image_.data[idx] == 0) {
        continue;
      }

      int left_label = (xx > 0) ? labels_[idx - 1] : 0;
      int top_label = (yy > 0) ? labels_[idx - width_] : 0;

      if (left_label == 0 && top_label == 0) {
        current_label_++;
        labels_[idx] = current_label_;
        parent_[current_label_] = current_label_;
      } else if (left_label != 0 && top_label == 0) {
        labels_[idx] = left_label;
      } else if (left_label == 0 && top_label != 0) {
        labels_[idx] = top_label;
      } else if (left_label == top_label) {
        labels_[idx] = left_label;
      } else {
        int min_label = std::min(left_label, top_label);
        int max_label = std::max(left_label, top_label);
        labels_[idx] = min_label;
        UnionLabels(min_label, max_label);
      }
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageSEQ::SecondPass() {
  std::unordered_map<int, int> new_labels;
  int next_label = 1;

  for (int ii = 0; ii < static_cast<int>(labels_.size()); ++ii) {
    if (labels_[ii] != 0) {
      int root = FindRoot(labels_[ii]);

      if (!new_labels.contains(root)) {
        new_labels[root] = next_label++;
      }

      labels_[ii] = new_labels[root];
    }
  }

  current_label_ = next_label - 1;
}

bool IvanovaPMarkingComponentsOnBinaryImageSEQ::RunImpl() {
  if (width_ <= 0 || height_ <= 0) {
    return false;
  }

  FirstPass();

  if (current_label_ > 0) {
    SecondPass();
  }

  return true;
}

bool IvanovaPMarkingComponentsOnBinaryImageSEQ::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();

  output.push_back(width_);
  output.push_back(height_);
  output.push_back(current_label_);

  for (int jj = 0; jj < static_cast<int>(labels_.size()); ++jj) {
    output.push_back(labels_[jj]);
  }

  return true;
}

}  // namespace ivanova_p_marking_components_on_binary_image
