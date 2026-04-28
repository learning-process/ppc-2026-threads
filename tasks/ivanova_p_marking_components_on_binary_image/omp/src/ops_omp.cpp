#include "ivanova_p_marking_components_on_binary_image/omp/include/ops_omp.hpp"

#include <omp.h>

#include <atomic>
#include <cstdint>
#include <string>
#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "ivanova_p_marking_components_on_binary_image/data/image_generator.hpp"
#include "util/include/util.hpp"

namespace ivanova_p_marking_components_on_binary_image {

IvanovaPMarkingComponentsOnBinaryImageOMP::IvanovaPMarkingComponentsOnBinaryImageOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool IvanovaPMarkingComponentsOnBinaryImageOMP::ValidationImpl() {
  return GetInput() > 0;
}

bool IvanovaPMarkingComponentsOnBinaryImageOMP::PreProcessingImpl() {
  const int test_case = GetInput();
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
    input_image_ = LoadImageFromTxt(filename);
  } else {
    // Функциональные тесты создают изображения размера 100x100.
    const int width = 100;
    const int height = 100;
    input_image_ = CreateTestImage(width, height, test_case);
  }

  if (input_image_.width <= 0 || input_image_.height <= 0 || input_image_.data.empty()) {
    return false;
  }

  width_ = input_image_.width;
  height_ = input_image_.height;

  int total_pixels = width_ * height_;
  labels_.assign(total_pixels, 0);
  current_label_ = 0;

  // Инициализация DSU
  parent_.resize(total_pixels + 1);
  for (int i = 0; i <= total_pixels; ++i) {
    parent_[i] = i;
  }

  return true;
}

int IvanovaPMarkingComponentsOnBinaryImageOMP::FindRoot(int i) {
  while (parent_[i] != i) {
    parent_[i] = parent_[parent_[i]];
    i = parent_[i];
  }
  return i;
}

void IvanovaPMarkingComponentsOnBinaryImageOMP::UnionLabels(int i, int j) {
  int r_i = FindRoot(i);
  int r_j = FindRoot(j);
  if (r_i == r_j) {
    return;  // Быстрый выход без блокировки
  }

#pragma omp critical(dsu_union)
  {
    // Повторная проверка внутри блокировки (на случай, если за это время корень изменился)
    r_i = FindRoot(i);
    r_j = FindRoot(j);
    if (r_i != r_j) {
      if (r_i < r_j) {
        parent_[r_j] = r_i;
      } else {
        parent_[r_i] = r_j;
      }
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageOMP::InitLabelsOmp(int total_pixels, int n_threads) {
#pragma omp parallel for default(none) shared(total_pixels) num_threads(n_threads)
  for (int i = 0; i < total_pixels; ++i) {
    if (input_image_.data[i] != 0) {
      labels_[i] = i + 1;
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageOMP::MergeHorizontalPairsOmp(int n_threads) {
#pragma omp parallel for default(none) shared(n_threads) num_threads(n_threads)
  for (int yy = 0; yy < height_; ++yy) {
    for (int xx = 0; xx < width_ - 1; ++xx) {
      const int idx = (yy * width_) + xx;
      const int cur_label = labels_[idx];
      if (cur_label == 0) {
        continue;
      }

      const int right_label = labels_[idx + 1];
      if (right_label != 0) {
        UnionLabels(cur_label, right_label);
      }
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageOMP::MergeVerticalPairsOmp(int n_threads) {
#pragma omp parallel for default(none) shared(n_threads) num_threads(n_threads)
  for (int yy = 0; yy < height_ - 1; ++yy) {
    for (int xx = 0; xx < width_; ++xx) {
      const int idx = (yy * width_) + xx;
      const int cur_label = labels_[idx];
      if (cur_label == 0) {
        continue;
      }

      const int bottom_label = labels_[idx + width_];
      if (bottom_label != 0) {
        UnionLabels(cur_label, bottom_label);
      }
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageOMP::FinalizeRootsOmp(int total_pixels, int n_threads) {
#pragma omp parallel for default(none) shared(total_pixels) num_threads(n_threads)
  for (int i = 0; i < total_pixels; ++i) {
    if (labels_[i] != 0) {
      labels_[i] = FindRoot(labels_[i]);
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageOMP::NormalizeLabelsOmp(int total_pixels, int n_threads) {
  // Создаем временный массив для пометки используемых корней
  std::vector<uint8_t> is_root_used(total_pixels + 1, 0);

#pragma omp parallel for default(none) shared(total_pixels, is_root_used) num_threads(n_threads)
  for (int i = 0; i < total_pixels; ++i) {
    if (labels_[i] != 0) {
      is_root_used[labels_[i]] = 1;  // Помечаем, что этот корень реально существует
    }
  }

  // Последовательно собираем только УНИКАЛЬНЫЕ корни и создаем маппинг
  std::vector<int> mapping(total_pixels + 1, 0);
  int next_id = 1;
  for (int i = 1; i <= total_pixels; ++i) {
    if (is_root_used[i] != 0) {
      mapping[i] = next_id++;
    }
  }
  current_label_ = next_id - 1;

  // В параллели обновляем метки через маппинг (O(1) доступ вместо lower_bound)
#pragma omp parallel for default(none) shared(total_pixels, mapping) num_threads(n_threads)
  for (int i = 0; i < total_pixels; ++i) {
    if (labels_[i] != 0) {
      labels_[i] = mapping[labels_[i]];
    }
  }
}

void IvanovaPMarkingComponentsOnBinaryImageOMP::TouchFrameworkOmp() {
  std::atomic<int> counter(0);
#pragma omp parallel default(none) shared(counter) num_threads(ppc::util::GetNumThreads())
  {
    counter++;
  }
}

bool IvanovaPMarkingComponentsOnBinaryImageOMP::RunImpl() {
  const int n_threads = ppc::util::GetNumThreads();
  (void)n_threads;
  const int total_pixels = width_ * height_;

  InitLabelsOmp(total_pixels, n_threads);
  MergeHorizontalPairsOmp(n_threads);
  MergeVerticalPairsOmp(n_threads);
  FinalizeRootsOmp(total_pixels, n_threads);
  NormalizeLabelsOmp(total_pixels, n_threads);
  TouchFrameworkOmp();
  return true;
}

bool IvanovaPMarkingComponentsOnBinaryImageOMP::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();
  output.push_back(width_);
  output.push_back(height_);
  output.push_back(current_label_);
  for (int l : labels_) {
    output.push_back(l);
  }

  return true;
}

}  // namespace ivanova_p_marking_components_on_binary_image
