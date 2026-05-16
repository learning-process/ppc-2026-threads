#include "marin_l_mark_components/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "marin_l_mark_components/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"

namespace marin_l_mark_components {

namespace {

constexpr std::uint64_t kMaxPixels = 100000000ULL;
constexpr int kMinRowsPerStripe = 64;

struct StripeSetup {
  int num_stripes = 1;
  int total_max_labels = 1;
  std::vector<int> stripe_bounds;
  std::vector<int> stripe_base_label;
};

std::vector<int> BuildCounts(int total_items, int proc_count) {
  std::vector<int> counts(static_cast<std::size_t>(proc_count), 0);
  for (int rank = 0; rank < proc_count; ++rank) {
    counts[static_cast<std::size_t>(rank)] =
        (((rank + 1) * total_items) / proc_count) - ((rank * total_items) / proc_count);
  }
  return counts;
}

std::vector<int> BuildDisplacements(const std::vector<int> &counts) {
  std::vector<int> displs(counts.size(), 0);
  for (std::size_t idx = 1; idx < counts.size(); ++idx) {
    displs[idx] = displs[idx - 1] + counts[idx - 1];
  }
  return displs;
}

std::vector<int> BuildPixelCounts(const std::vector<int> &row_counts, int width) {
  std::vector<int> pixel_counts(row_counts.size(), 0);
  for (std::size_t rank = 0; rank < row_counts.size(); ++rank) {
    pixel_counts[rank] = row_counts[rank] * width;
  }
  return pixel_counts;
}

std::vector<std::uint8_t> FlattenImage(const Image &img) {
  if (img.empty() || img.front().empty()) {
    return {};
  }

  const int height = static_cast<int>(img.size());
  const int width = static_cast<int>(img.front().size());
  std::vector<std::uint8_t> flat(static_cast<std::size_t>(height) * static_cast<std::size_t>(width), 0);

  for (int row = 0; row < height; ++row) {
    const std::size_t row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width);
    for (int col = 0; col < width; ++col) {
      flat[row_offset + static_cast<std::size_t>(col)] =
          static_cast<std::uint8_t>(img[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]);
    }
  }

  return flat;
}

int FindRoot(std::vector<int> &parent, int x) {
  int root = x;
  while (parent[static_cast<std::size_t>(root)] != root) {
    root = parent[static_cast<std::size_t>(root)];
  }

  int current = x;
  while (current != root) {
    const int next = parent[static_cast<std::size_t>(current)];
    parent[static_cast<std::size_t>(current)] = root;
    current = next;
  }

  return root;
}

void UnionLabels(std::vector<int> &parent, int a, int b) {
  const int root_a = FindRoot(parent, a);
  const int root_b = FindRoot(parent, b);
  if (root_a == root_b) {
    return;
  }

  if (root_a < root_b) {
    parent[static_cast<std::size_t>(root_b)] = root_a;
  } else {
    parent[static_cast<std::size_t>(root_a)] = root_b;
  }
}

StripeSetup BuildStripeSetup(int height, int width, int concurrency) {
  StripeSetup setup;
  if (height <= 0 || width <= 0) {
    return setup;
  }

  setup.num_stripes = std::min(height, concurrency * 2);
  if (setup.num_stripes > 0 && height / setup.num_stripes < kMinRowsPerStripe) {
    setup.num_stripes = std::max(1, height / kMinRowsPerStripe);
  }
  setup.num_stripes = std::max(1, setup.num_stripes);

  setup.stripe_bounds.assign(static_cast<std::size_t>(setup.num_stripes) + 1ULL, 0);
  for (int stripe = 0; stripe <= setup.num_stripes; ++stripe) {
    setup.stripe_bounds[static_cast<std::size_t>(stripe)] = (stripe * height) / setup.num_stripes;
  }

  setup.stripe_base_label.assign(static_cast<std::size_t>(setup.num_stripes), 0);
  for (int stripe = 0; stripe < setup.num_stripes; ++stripe) {
    setup.stripe_base_label[static_cast<std::size_t>(stripe)] = setup.total_max_labels;
    const int stripe_height = setup.stripe_bounds[static_cast<std::size_t>(stripe) + 1ULL] -
                              setup.stripe_bounds[static_cast<std::size_t>(stripe)];
    setup.total_max_labels += ((stripe_height * width) / 2) + 1;
  }

  return setup;
}

void InitializeParents(std::vector<int> &parent, int total_max_labels) {
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, total_max_labels),
                            [&](const oneapi::tbb::blocked_range<int> &range) {
    for (int label = range.begin(); label != range.end(); ++label) {
      parent[static_cast<std::size_t>(label)] = label;
    }
  });
}

void AssignPixelLabel(std::vector<int> &labels_flat, std::vector<int> &parent, std::size_t idx, int left_label,
                      int top_label, int &next_label) {
  if (left_label == 0 && top_label == 0) {
    labels_flat[idx] = next_label++;
    return;
  }
  if (left_label == 0) {
    labels_flat[idx] = top_label;
    return;
  }
  if (top_label == 0) {
    labels_flat[idx] = left_label;
    return;
  }

  const int min_label = std::min(left_label, top_label);
  labels_flat[idx] = min_label;
  if (left_label != top_label) {
    UnionLabels(parent, left_label, top_label);
  }
}

void LabelStripe(const std::vector<std::uint8_t> &binary_flat, std::vector<int> &labels_flat, std::vector<int> &parent,
                 const std::vector<int> &stripe_bounds, const std::vector<int> &stripe_base_label,
                 std::vector<int> &stripe_label_end, int width, int stripe) {
  const int start_row = stripe_bounds[static_cast<std::size_t>(stripe)];
  const int end_row = stripe_bounds[static_cast<std::size_t>(stripe) + 1ULL];
  int next_label = stripe_base_label[static_cast<std::size_t>(stripe)];

  for (int row = start_row; row < end_row; ++row) {
    const std::size_t row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width);
    const bool has_top = row > start_row;
    const std::size_t top_row_offset = has_top ? row_offset - static_cast<std::size_t>(width) : 0;

    for (int col = 0; col < width; ++col) {
      const std::size_t idx = row_offset + static_cast<std::size_t>(col);
      if (binary_flat[idx] == 0U) {
        continue;
      }

      const int left_label = (col > 0) ? labels_flat[idx - 1ULL] : 0;
      const int top_label = has_top ? labels_flat[top_row_offset + static_cast<std::size_t>(col)] : 0;
      AssignPixelLabel(labels_flat, parent, idx, left_label, top_label, next_label);
    }
  }

  stripe_label_end[static_cast<std::size_t>(stripe)] = next_label;
}

void MergeLocalStripeBorders(const std::vector<std::uint8_t> &binary_flat, const std::vector<int> &stripe_bounds,
                             std::vector<int> &labels_flat, std::vector<int> &parent, int width, int num_stripes) {
  for (int stripe = 0; stripe < num_stripes - 1; ++stripe) {
    const int boundary_row = stripe_bounds[static_cast<std::size_t>(stripe) + 1ULL];
    const std::size_t top_row_offset = static_cast<std::size_t>(boundary_row - 1) * static_cast<std::size_t>(width);
    const std::size_t bottom_row_offset = static_cast<std::size_t>(boundary_row) * static_cast<std::size_t>(width);

    for (int col = 0; col < width; ++col) {
      const std::size_t top_idx = top_row_offset + static_cast<std::size_t>(col);
      const std::size_t bottom_idx = bottom_row_offset + static_cast<std::size_t>(col);
      if (binary_flat[top_idx] == 1U && binary_flat[bottom_idx] == 1U) {
        const int top_label = labels_flat[top_idx];
        const int bottom_label = labels_flat[bottom_idx];
        if (top_label > 0 && bottom_label > 0 && top_label != bottom_label) {
          UnionLabels(parent, top_label, bottom_label);
        }
      }
    }
  }
}

int CompactLocalLabels(std::vector<int> &labels_flat, std::vector<int> &parent,
                       const std::vector<int> &stripe_base_label, const std::vector<int> &stripe_label_end,
                       int total_max_labels) {
  for (int label = 1; label < total_max_labels; ++label) {
    parent[static_cast<std::size_t>(label)] = FindRoot(parent, label);
  }

  std::vector<int> compacted(static_cast<std::size_t>(total_max_labels), 0);
  int next_compact_id = 1;

  for (std::size_t stripe = 0; stripe < stripe_base_label.size(); ++stripe) {
    for (int label = stripe_base_label[stripe]; label < stripe_label_end[stripe]; ++label) {
      const int root = parent[static_cast<std::size_t>(label)];
      if (compacted[static_cast<std::size_t>(root)] == 0) {
        compacted[static_cast<std::size_t>(root)] = next_compact_id++;
      }
      compacted[static_cast<std::size_t>(label)] = compacted[static_cast<std::size_t>(root)];
    }
  }

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, labels_flat.size()),
                            [&](const oneapi::tbb::blocked_range<std::size_t> &range) {
    for (std::size_t idx = range.begin(); idx != range.end(); ++idx) {
      const int label = labels_flat[idx];
      if (label != 0) {
        labels_flat[idx] = compacted[static_cast<std::size_t>(label)];
      }
    }
  });

  return next_compact_id - 1;
}

void ApplyRankOffset(std::vector<int> &labels_flat, int label_offset) {
  if (label_offset == 0) {
    return;
  }

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, labels_flat.size()),
                            [&](const oneapi::tbb::blocked_range<std::size_t> &range) {
    for (std::size_t idx = range.begin(); idx != range.end(); ++idx) {
      if (labels_flat[idx] != 0) {
        labels_flat[idx] += label_offset;
      }
    }
  });
}

std::size_t FindNextNonEmptyRank(const std::vector<int> &row_counts, std::size_t rank) {
  std::size_t next_rank = rank + 1;
  while (next_rank < row_counts.size() && row_counts[next_rank] == 0) {
    ++next_rank;
  }
  return next_rank;
}

bool HasAdjacentBoundary(const std::vector<int> &row_counts, const std::vector<int> &row_displs, std::size_t rank,
                         std::size_t next_rank) {
  if (next_rank >= row_counts.size()) {
    return false;
  }

  const int boundary_row = row_displs[rank] + row_counts[rank];
  return row_displs[next_rank] == boundary_row;
}

void MergeBoundaryRow(const std::vector<std::uint8_t> &global_binary_flat, std::vector<int> &global_labels_flat,
                      std::vector<int> &parent, int boundary_row, int width) {
  const std::size_t top_row_offset = static_cast<std::size_t>(boundary_row - 1) * static_cast<std::size_t>(width);
  const std::size_t bottom_row_offset = static_cast<std::size_t>(boundary_row) * static_cast<std::size_t>(width);

  for (int col = 0; col < width; ++col) {
    const std::size_t top_idx = top_row_offset + static_cast<std::size_t>(col);
    const std::size_t bottom_idx = bottom_row_offset + static_cast<std::size_t>(col);
    if (global_binary_flat[top_idx] != 1U || global_binary_flat[bottom_idx] != 1U) {
      continue;
    }

    const int top_label = global_labels_flat[top_idx];
    const int bottom_label = global_labels_flat[bottom_idx];
    if (top_label > 0 && bottom_label > 0 && top_label != bottom_label) {
      UnionLabels(parent, top_label, bottom_label);
    }
  }
}

void MergeRankBorders(const std::vector<std::uint8_t> &global_binary_flat, std::vector<int> &global_labels_flat,
                      std::vector<int> &parent, const std::vector<int> &row_counts, const std::vector<int> &row_displs,
                      int width) {
  for (std::size_t rank = 0; rank < row_counts.size(); ++rank) {
    if (row_counts[rank] == 0) {
      continue;
    }

    const std::size_t next_rank = FindNextNonEmptyRank(row_counts, rank);
    if (!HasAdjacentBoundary(row_counts, row_displs, rank, next_rank)) {
      if (next_rank >= row_counts.size()) {
        break;
      }
      continue;
    }

    const int boundary_row = row_displs[rank] + row_counts[rank];
    if (boundary_row <= 0) {
      break;
    }
    MergeBoundaryRow(global_binary_flat, global_labels_flat, parent, boundary_row, width);
  }
}

void CompactGlobalLabels(std::vector<int> &global_labels_flat, std::vector<int> &parent, int total_labels) {
  for (int label = 1; label <= total_labels; ++label) {
    parent[static_cast<std::size_t>(label)] = FindRoot(parent, label);
  }

  std::vector<int> compacted(static_cast<std::size_t>(total_labels) + 1ULL, 0);
  int next_compact_id = 1;

  for (int &label : global_labels_flat) {
    if (label == 0) {
      continue;
    }

    const int root = parent[static_cast<std::size_t>(label)];
    if (compacted[static_cast<std::size_t>(root)] == 0) {
      compacted[static_cast<std::size_t>(root)] = next_compact_id++;
    }
    label = compacted[static_cast<std::size_t>(root)];
  }
}

}  // namespace

MarinLMarkComponentsALL::MarinLMarkComponentsALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool MarinLMarkComponentsALL::IsBinary(const Image &img) {
  for (const auto &row : img) {
    for (int pixel : row) {
      if (pixel != 0 && pixel != 1) {
        return false;
      }
    }
  }
  return true;
}

bool MarinLMarkComponentsALL::ValidationImpl() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);

  const auto &img = GetInput().binary;
  if (img.empty() || img.front().empty()) {
    return false;
  }

  const std::size_t width = img.front().size();
  for (const auto &row : img) {
    if (row.size() != width) {
      return false;
    }
  }

  return IsBinary(img);
}

bool MarinLMarkComponentsALL::PreProcessingImpl() {
  const auto &img = GetInput().binary;
  height_ = static_cast<int>(img.size());
  width_ = static_cast<int>(img.front().size());

  if (height_ <= 0 || width_ <= 0) {
    return false;
  }

  const std::uint64_t total_pixels = static_cast<std::uint64_t>(height_) * static_cast<std::uint64_t>(width_);
  if (total_pixels > kMaxPixels) {
    return false;
  }

  local_binary_flat_.clear();
  global_labels_flat_.assign(static_cast<std::size_t>(total_pixels), 0);
  labels_out_.clear();
  return true;
}

bool MarinLMarkComponentsALL::RunImpl() {
  const std::vector<int> row_counts = BuildCounts(height_, world_size_);
  const std::vector<int> row_displs = BuildDisplacements(row_counts);
  const std::vector<int> pixel_counts = BuildPixelCounts(row_counts, width_);
  const std::vector<int> pixel_displs = BuildDisplacements(pixel_counts);

  std::vector<std::uint8_t> global_binary_flat;
  if (rank_ == 0) {
    global_binary_flat = FlattenImage(GetInput().binary);
  }

  const int local_height = row_counts[static_cast<std::size_t>(rank_)];
  const int local_pixels = pixel_counts[static_cast<std::size_t>(rank_)];
  local_binary_flat_.assign(static_cast<std::size_t>(local_pixels), 0);
  std::vector<int> local_labels_flat(static_cast<std::size_t>(local_pixels), 0);

  MPI_Scatterv(global_binary_flat.data(), pixel_counts.data(), pixel_displs.data(), MPI_UINT8_T,
               local_binary_flat_.data(), local_pixels, MPI_UINT8_T, 0, MPI_COMM_WORLD);

  const int thread_count = std::max(1, ppc::util::GetNumThreads());
  const oneapi::tbb::global_control control(oneapi::tbb::global_control::max_allowed_parallelism,
                                            static_cast<std::size_t>(thread_count));

  int local_component_count = 0;
  if (local_height > 0) {
    const StripeSetup setup = BuildStripeSetup(local_height, width_, thread_count);
    std::vector<int> parent(static_cast<std::size_t>(setup.total_max_labels), 0);
    InitializeParents(parent, setup.total_max_labels);

    std::vector<int> stripe_label_end(static_cast<std::size_t>(setup.num_stripes), 0);
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, setup.num_stripes),
                              [&](const oneapi::tbb::blocked_range<int> &range) {
      for (int stripe = range.begin(); stripe != range.end(); ++stripe) {
        LabelStripe(local_binary_flat_, local_labels_flat, parent, setup.stripe_bounds, setup.stripe_base_label,
                    stripe_label_end, width_, stripe);
      }
    });

    MergeLocalStripeBorders(local_binary_flat_, setup.stripe_bounds, local_labels_flat, parent, width_,
                            setup.num_stripes);
    local_component_count = CompactLocalLabels(local_labels_flat, parent, setup.stripe_base_label, stripe_label_end,
                                               setup.total_max_labels);
  }

  std::vector<int> local_component_counts(static_cast<std::size_t>(world_size_), 0);
  MPI_Allgather(&local_component_count, 1, MPI_INT, local_component_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> component_offsets(static_cast<std::size_t>(world_size_), 0);
  int total_component_count = 0;
  for (int rank = 0; rank < world_size_; ++rank) {
    component_offsets[static_cast<std::size_t>(rank)] = total_component_count;
    total_component_count += local_component_counts[static_cast<std::size_t>(rank)];
  }

  ApplyRankOffset(local_labels_flat, component_offsets[static_cast<std::size_t>(rank_)]);

  MPI_Gatherv(local_labels_flat.data(), local_pixels, MPI_INT, global_labels_flat_.data(), pixel_counts.data(),
              pixel_displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank_ == 0) {
    std::vector<int> global_parent(static_cast<std::size_t>(total_component_count) + 1ULL, 0);
    std::ranges::generate(global_parent, [value = 0]() mutable { return value++; });
    MergeRankBorders(global_binary_flat, global_labels_flat_, global_parent, row_counts, row_displs, width_);
    CompactGlobalLabels(global_labels_flat_, global_parent, total_component_count);
  }

  MPI_Bcast(global_labels_flat_.data(), static_cast<int>(global_labels_flat_.size()), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool MarinLMarkComponentsALL::PostProcessingImpl() {
  ConvertLabelsToOutput();
  OutType out;
  out.labels = labels_out_;
  GetOutput() = out;
  return true;
}

void MarinLMarkComponentsALL::ConvertLabelsToOutput() {
  labels_out_.resize(static_cast<std::size_t>(height_));

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, height_),
                            [&](const oneapi::tbb::blocked_range<int> &range) {
    for (int row = range.begin(); row != range.end(); ++row) {
      labels_out_[static_cast<std::size_t>(row)].resize(static_cast<std::size_t>(width_));
      const std::size_t row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width_);
      for (int col = 0; col < width_; ++col) {
        labels_out_[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] =
            global_labels_flat_[row_offset + static_cast<std::size_t>(col)];
      }
    }
  });
}

}  // namespace marin_l_mark_components
