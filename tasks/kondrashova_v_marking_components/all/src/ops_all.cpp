#include "kondrashova_v_marking_components/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "kondrashova_v_marking_components/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kondrashova_v_marking_components {

KondrashovaVTaskALL::KondrashovaVTaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool KondrashovaVTaskALL::ValidationImpl() {
  return true;
}

bool KondrashovaVTaskALL::PreProcessingImpl() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);

  if (rank_ == 0) {
    const auto &in = GetInput();
    width_ = in.width;
    height_ = in.height;
    image_ = in.data;
  }

  MPI_Bcast(&width_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&height_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int has_valid_input = 0;
  if (rank_ == 0) {
    has_valid_input = (width_ > 0 && height_ > 0 && static_cast<int>(image_.size()) == (width_ * height_)) ? 1 : 0;
  }
  MPI_Bcast(&has_valid_input, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (has_valid_input == 0) {
    image_.clear();
    labels_1d_.clear();
    GetOutput().count = 0;
    GetOutput().labels.clear();
    return true;
  }

  if (rank_ != 0) {
    image_.resize(static_cast<size_t>(width_) * static_cast<size_t>(height_));
  }

  MPI_Bcast(image_.data(), width_ * height_, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  labels_1d_.assign(static_cast<size_t>(width_) * static_cast<size_t>(height_), 0);

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

void MergeVertical(int height, int width, const std::vector<int> &local_labels, std::vector<int> &parent,
                   std::vector<int> &rnk) {
  for (int ii = 1; ii < height; ++ii) {
    for (int jj = 0; jj < width; ++jj) {
      auto idx = (static_cast<size_t>(ii) * static_cast<size_t>(width)) + static_cast<size_t>(jj);
      auto tidx = (static_cast<size_t>(ii - 1) * static_cast<size_t>(width)) + static_cast<size_t>(jj);
      if (local_labels[idx] != 0 && local_labels[tidx] != 0 && local_labels[idx] != local_labels[tidx]) {
        Unite(parent, rnk, local_labels[idx], local_labels[tidx]);
      }
    }
  }
}

void MergeBoundaries(int row_start, int row_end, int width, int num_threads, const std::vector<int> &local_labels,
                     std::vector<int> &parent, std::vector<int> &rnk) {
  const int rows = row_end - row_start;
  for (int tid = 1; tid < num_threads; ++tid) {
    const int boundary_row = row_start + ((tid * rows) / num_threads);
    if (boundary_row <= row_start || boundary_row >= row_end) {
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

void FillRowLabels(int width, int row, const std::vector<int> &local_labels, std::vector<int> &row_labels) {
  for (int jj = 0; jj < width; ++jj) {
    row_labels[static_cast<size_t>(jj)] =
        local_labels[(static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(jj)];
  }
}

void MergeBoundaryLabels(int width, const std::vector<int> &send_row, const std::vector<int> &recv_row,
                         std::vector<int> &parent, std::vector<int> &rnk) {
  for (int jj = 0; jj < width; ++jj) {
    const int label_cur = send_row[static_cast<size_t>(jj)];
    const int label_neighbor = recv_row[static_cast<size_t>(jj)];
    if (label_cur != 0 && label_neighbor != 0 && label_cur != label_neighbor) {
      Unite(parent, rnk, label_cur, label_neighbor);
    }
  }
}

void ExchangeAndMergeRow(int width, int neighbor_rank, int send_tag, int recv_tag, bool has_rows, int row_index,
                         const std::vector<int> &local_labels, std::vector<int> &parent, std::vector<int> &rnk) {
  std::vector<int> send_row(static_cast<size_t>(width), 0);
  std::vector<int> recv_row(static_cast<size_t>(width));

  if (has_rows) {
    FillRowLabels(width, row_index, local_labels, send_row);
  }

  MPI_Sendrecv(send_row.data(), width, MPI_INT, neighbor_rank, send_tag, recv_row.data(), width, MPI_INT, neighbor_rank,
               recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  MergeBoundaryLabels(width, send_row, recv_row, parent, rnk);
}

void MergeMPIBoundaries(int width, int rank, int world_size, int local_row_start, int local_row_end,
                        const std::vector<int> &local_labels, std::vector<int> &parent, std::vector<int> &rnk) {
  const bool has_rows = local_row_start < local_row_end;
  const int first_row = local_row_start;
  const int last_row = has_rows ? (local_row_end - 1) : local_row_start;

  if (rank > 0) {
    ExchangeAndMergeRow(width, rank - 1, 0, 1, has_rows, first_row, local_labels, parent, rnk);
  }

  if (rank < world_size - 1) {
    ExchangeAndMergeRow(width, rank + 1, 1, 0, has_rows, last_row, local_labels, parent, rnk);
  }
}

}  // namespace

bool KondrashovaVTaskALL::RunImpl() {
  if (width_ <= 0 || height_ <= 0 || image_.empty()) {
    GetOutput().count = 0;
    return true;
  }

  const int total = width_ * height_;
  const int num_threads = ppc::util::GetNumThreads();

  const int mpi_row_start = (rank_ * height_) / world_size_;
  const int mpi_row_end = ((rank_ + 1) * height_) / world_size_;
  const int mpi_rows = mpi_row_end - mpi_row_start;

  const int max_rows_per_rank = (height_ + world_size_ - 1) / world_size_;
  const int max_labels_per_thread = (max_rows_per_rank * width_) + 1;

  std::vector<int> local_labels(static_cast<size_t>(total), 0);

  const int width = width_;
  const std::vector<uint8_t> image = image_;

#pragma omp parallel num_threads(num_threads) default(none) shared(local_labels, width, image, num_threads) \
    firstprivate(mpi_row_start, mpi_rows, max_labels_per_thread)
  {
    const int tid = omp_get_thread_num();
    const int omp_row_start = mpi_row_start + ((tid * mpi_rows) / num_threads);
    const int omp_row_end = mpi_row_start + (((tid + 1) * mpi_rows) / num_threads);
    const int label_offset = (rank_ * num_threads * max_labels_per_thread) + (tid * max_labels_per_thread);
    ScanStripe(omp_row_start, omp_row_end, width, label_offset, image, local_labels);
  }

  const int global_max_labels = (world_size_ * num_threads * max_labels_per_thread) + 1;
  std::vector<int> parent(static_cast<size_t>(global_max_labels));
  std::vector<int> rnk(static_cast<size_t>(global_max_labels), 0);
  for (int ii = 0; ii < global_max_labels; ++ii) {
    parent[static_cast<size_t>(ii)] = ii;
  }

  MergeHorizontal(height_, width_, local_labels, parent, rnk);
  MergeVertical(height_, width_, local_labels, parent, rnk);
  MergeBoundaries(mpi_row_start, mpi_row_end, width_, num_threads, local_labels, parent, rnk);

  MergeMPIBoundaries(width_, rank_, world_size_, mpi_row_start, mpi_row_end, local_labels, parent, rnk);

  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<int> all_local_labels(static_cast<size_t>(total), 0);
  MPI_Allreduce(local_labels.data(), all_local_labels.data(), total, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (rank_ == 0) {
    std::vector<int> global_parent(static_cast<size_t>(global_max_labels));
    std::vector<int> global_rnk(static_cast<size_t>(global_max_labels), 0);
    for (int ii = 0; ii < global_max_labels; ++ii) {
      global_parent[static_cast<size_t>(ii)] = ii;
    }

    MergeHorizontal(height_, width_, all_local_labels, global_parent, global_rnk);
    MergeVertical(height_, width_, all_local_labels, global_parent, global_rnk);

    std::vector<int> relabel_map(static_cast<size_t>(global_max_labels), 0);
    GetOutput().count = Relabel(total, all_local_labels, global_parent, relabel_map, labels_1d_);
  }

  int global_count = GetOutput().count;
  MPI_Bcast(&global_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  GetOutput().count = global_count;
  if (total > 0) {
    if (labels_1d_.size() != static_cast<size_t>(total)) {
      labels_1d_.assign(static_cast<size_t>(total), 0);
    }
    MPI_Bcast(labels_1d_.data(), total, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return true;
}

bool KondrashovaVTaskALL::PostProcessingImpl() {
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
