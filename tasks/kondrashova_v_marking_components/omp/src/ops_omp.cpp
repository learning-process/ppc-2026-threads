#include "kondrashova_v_marking_components/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "kondrashova_v_marking_components/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kondrashova_v_marking_components {

KondrashovaVTaskOMP::KondrashovaVTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool KondrashovaVTaskOMP::ValidationImpl() {
  return true;
}

bool KondrashovaVTaskOMP::PreProcessingImpl() {
  const auto &in = GetInput();

  image_ = in.data;
  width_ = in.width;
  height_ = in.height;

  if (width_ > 0 && height_ > 0 && static_cast<int>(image_.size()) == width_ * height_) {
    labels_1d_.assign(static_cast<size_t>(width_) * static_cast<size_t>(height_), 0);
  } else {
    labels_1d_.clear();
  }

  GetOutput().count = 0;
  GetOutput().labels.clear();
  return true;
}

namespace {

int Find(std::vector<int> &parent, int xx) {
  while (parent[static_cast<size_t>(xx)] != xx) {
    parent[static_cast<size_t>(xx)] = parent[static_cast<size_t>(parent[static_cast<size_t>(xx)])];
    xx = parent[static_cast<size_t>(xx)];
  }
  return xx;
}

void Unite(std::vector<int> &parent, std::vector<int> &rnk, int aa, int bb) {
  aa = Find(parent, aa);
  bb = Find(parent, bb);
  if (aa == bb) {
    return;
  }
  if (rnk[static_cast<size_t>(aa)] < rnk[static_cast<size_t>(bb)]) {
    std::swap(aa, bb);
  }
  parent[static_cast<size_t>(bb)] = aa;
  if (rnk[static_cast<size_t>(aa)] == rnk[static_cast<size_t>(bb)]) {
    rnk[static_cast<size_t>(aa)]++;
  }
}

int GetNeighborLabel(int ii, int jj, int di, int dj, int row_start, int row_end, int width,
                     const std::vector<uint8_t> &image, const std::vector<int> &local_labels) {
  int ni = ii + di;
  int nj = jj + dj;
  if (ni < row_start || ni >= row_end || nj < 0 || nj >= width) {
    return 0;
  }
  auto nidx = (static_cast<size_t>(ni) * static_cast<size_t>(width)) + static_cast<size_t>(nj);
  if (image[nidx] == 0) {
    return local_labels[nidx];
  }
  return 0;
}

void ScanStripe(int row_start, int row_end, int width, int label_offset, const std::vector<uint8_t> &image,
                std::vector<int> &local_labels) {
  int current_label = label_offset;
  for (int ii = row_start; ii < row_end; ++ii) {
    for (int jj = 0; jj < width; ++jj) {
      auto idx = (static_cast<size_t>(ii) * static_cast<size_t>(width)) + static_cast<size_t>(jj);
      if (image[idx] != 0) {
        continue;
      }

      int left_label = GetNeighborLabel(ii, jj, 0, -1, row_start, row_end, width, image, local_labels);
      int top_label = GetNeighborLabel(ii, jj, -1, 0, row_start, row_end, width, image, local_labels);

      if (left_label == 0 && top_label == 0) {
        local_labels[idx] = ++current_label;
      } else if (left_label != 0 && top_label == 0) {
        local_labels[idx] = left_label;
      } else if (left_label == 0) {
        local_labels[idx] = top_label;
      } else {
        local_labels[idx] = std::min(left_label, top_label);
      }
    }
  }
}

void MergeHorizontal(int height, int width, const std::vector<int> &local_labels, std::vector<int> &parent,
                     std::vector<int> &rnk) {
  for (int ii = 0; ii < height; ++ii) {
    for (int jj = 1; jj < width; ++jj) {
      auto idx = (static_cast<size_t>(ii) * static_cast<size_t>(width)) + static_cast<size_t>(jj);
      auto lidx = (static_cast<size_t>(ii) * static_cast<size_t>(width)) + static_cast<size_t>(jj - 1);
      if (local_labels[idx] != 0 && local_labels[lidx] != 0 && local_labels[idx] != local_labels[lidx]) {
        Unite(parent, rnk, local_labels[idx], local_labels[lidx]);
      }
    }
  }
}

void MergeBoundaries(int height, int width, int num_threads, const std::vector<int> &local_labels,
                     std::vector<int> &parent, std::vector<int> &rnk) {
  for (int tid = 1; tid < num_threads; ++tid) {
    const int boundary_row = (tid * height) / num_threads;
    if (boundary_row >= height) {
      continue;
    }
    for (int jj = 0; jj < width; ++jj) {
      auto idx = (static_cast<size_t>(boundary_row) * static_cast<size_t>(width)) + static_cast<size_t>(jj);
      auto tidx = (static_cast<size_t>(boundary_row - 1) * static_cast<size_t>(width)) + static_cast<size_t>(jj);
      if (local_labels[idx] != 0 && local_labels[tidx] != 0 && local_labels[idx] != local_labels[tidx]) {
        Unite(parent, rnk, local_labels[idx], local_labels[tidx]);
      }
    }
  }
}

int Relabel(int total, const std::vector<int> &local_labels, std::vector<int> &parent, std::vector<int> &relabel_map,
            std::vector<int> &labels_1d) {
  int count = 0;
  for (int ii = 0; ii < total; ++ii) {
    auto idx = static_cast<size_t>(ii);
    if (local_labels[idx] == 0) {
      continue;
    }
    int root = Find(parent, local_labels[idx]);
    if (relabel_map[static_cast<size_t>(root)] == 0) {
      relabel_map[static_cast<size_t>(root)] = ++count;
    }
    labels_1d[idx] = relabel_map[static_cast<size_t>(root)];
  }
  return count;
}

}  // namespace

bool KondrashovaVTaskOMP::RunImpl() {
  if (width_ <= 0 || height_ <= 0 || image_.empty()) {
    GetOutput().count = 0;
    return true;
  }

  const int total = width_ * height_;
  const int num_threads = ppc::util::GetNumThreads();
  const int max_labels_per_thread = total + 1;
  const int max_total_labels = (num_threads * max_labels_per_thread) + 1;

  std::vector<int> local_labels(static_cast<size_t>(total), 0);

  // Копируем члены класса в локальные переменные для OpenMP (требование MSVC)
  const int width = width_;
  const int height = height_;
  const std::vector<uint8_t> image = image_;

#pragma omp parallel num_threads(num_threads) default(none) shared(local_labels, width, height, image, num_threads) \
    firstprivate(max_labels_per_thread)
  {
    const int tid = omp_get_thread_num();
    const int row_start = (tid * height) / num_threads;
    const int row_end = ((tid + 1) * height) / num_threads;
    const int label_offset = tid * max_labels_per_thread;
    ScanStripe(row_start, row_end, width, label_offset, image, local_labels);
  }

  std::vector<int> parent(static_cast<size_t>(max_total_labels));
  std::vector<int> rnk(static_cast<size_t>(max_total_labels), 0);
  for (int ii = 0; ii < max_total_labels; ++ii) {
    parent[static_cast<size_t>(ii)] = ii;
  }

  MergeHorizontal(height, width, local_labels, parent, rnk);
  MergeBoundaries(height, width, num_threads, local_labels, parent, rnk);

  std::vector<int> relabel_map(static_cast<size_t>(max_total_labels), 0);
  GetOutput().count = Relabel(total, local_labels, parent, relabel_map, labels_1d_);

  return true;
}

bool KondrashovaVTaskOMP::PostProcessingImpl() {
  if (width_ <= 0 || height_ <= 0) {
    GetOutput().labels.clear();
    return true;
  }

  GetOutput().labels.assign(height_, std::vector<int>(width_, 0));
  for (int ii = 0; ii < height_; ++ii) {
    for (int jj = 0; jj < width_; ++jj) {
      auto idx = (static_cast<size_t>(ii) * static_cast<size_t>(width_)) + static_cast<size_t>(jj);
      GetOutput().labels[ii][jj] = labels_1d_[idx];
    }
  }
  return true;
}

}  // namespace kondrashova_v_marking_components
