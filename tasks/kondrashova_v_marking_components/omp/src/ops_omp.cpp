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

  width_ = in.width;
  height_ = in.height;
  image_ = in.data;

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
  while (parent[xx] != xx) {
    parent[xx] = parent[parent[xx]];
    xx = parent[xx];
  }
  return xx;
}

void Unite(std::vector<int> &parent, std::vector<int> &rnk, int aa, int bb) {
  aa = Find(parent, aa);
  bb = Find(parent, bb);
  if (aa == bb) {
    return;
  }
  if (rnk[aa] < rnk[bb]) {
    std::swap(aa, bb);
  }
  parent[bb] = aa;
  if (rnk[aa] == rnk[bb]) {
    rnk[aa]++;
  }
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

#pragma omp parallel num_threads(num_threads)
  {
    const int tid = omp_get_thread_num();
    const int row_start = (tid * height_) / num_threads;
    const int row_end = ((tid + 1) * height_) / num_threads;
    const int label_offset = tid * max_labels_per_thread;
    int current_label = label_offset;

    for (int ii = row_start; ii < row_end; ++ii) {
      for (int jj = 0; jj < width_; ++jj) {
        auto idx = (static_cast<size_t>(ii) * static_cast<size_t>(width_)) + static_cast<size_t>(jj);
        if (image_[idx] != 0) {
          continue;
        }

        int left_label = 0;
        int top_label = 0;

        if (jj > 0) {
          auto lidx = (static_cast<size_t>(ii) * static_cast<size_t>(width_)) + static_cast<size_t>(jj - 1);
          if (image_[lidx] == 0) {
            left_label = local_labels[lidx];
          }
        }
        if (ii > row_start) {
          auto tidx = (static_cast<size_t>(ii - 1) * static_cast<size_t>(width_)) + static_cast<size_t>(jj);
          if (image_[tidx] == 0) {
            top_label = local_labels[tidx];
          }
        }

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

  std::vector<int> parent(static_cast<size_t>(max_total_labels));
  std::vector<int> rnk(static_cast<size_t>(max_total_labels), 0);
  for (int ii = 0; ii < max_total_labels; ++ii) {
    parent[static_cast<size_t>(ii)] = ii;
  }

  for (int ii = 0; ii < height_; ++ii) {
    for (int jj = 1; jj < width_; ++jj) {
      auto idx = (static_cast<size_t>(ii) * static_cast<size_t>(width_)) + static_cast<size_t>(jj);
      auto lidx = (static_cast<size_t>(ii) * static_cast<size_t>(width_)) + static_cast<size_t>(jj - 1);
      if (local_labels[idx] != 0 && local_labels[lidx] != 0 && local_labels[idx] != local_labels[lidx]) {
        Unite(parent, rnk, local_labels[idx], local_labels[lidx]);
      }
    }
  }

  for (int tid = 1; tid < num_threads; ++tid) {
    const int boundary_row = (tid * height_) / num_threads;
    if (boundary_row >= height_) {
      continue;
    }
    for (int jj = 0; jj < width_; ++jj) {
      auto idx = (static_cast<size_t>(boundary_row) * static_cast<size_t>(width_)) + static_cast<size_t>(jj);
      auto tidx = (static_cast<size_t>(boundary_row - 1) * static_cast<size_t>(width_)) + static_cast<size_t>(jj);
      if (local_labels[idx] != 0 && local_labels[tidx] != 0 && local_labels[idx] != local_labels[tidx]) {
        Unite(parent, rnk, local_labels[idx], local_labels[tidx]);
      }
    }
  }

  std::vector<int> relabel(static_cast<size_t>(max_total_labels), 0);
  int count = 0;
  for (int ii = 0; ii < total; ++ii) {
    auto idx = static_cast<size_t>(ii);
    if (local_labels[idx] == 0) {
      continue;
    }
    int root = Find(parent, local_labels[idx]);
    if (relabel[static_cast<size_t>(root)] == 0) {
      relabel[static_cast<size_t>(root)] = ++count;
    }
    labels_1d_[idx] = relabel[static_cast<size_t>(root)];
  }

  GetOutput().count = count;
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
