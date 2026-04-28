#include "marin_l_mark_components/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "marin_l_mark_components/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/task_arena.h"

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

int FindRoot(std::vector<int> &parent, int x) {
  int root = x;
  while (parent[root] != root) {
    root = parent[root];
  }

  int current = x;
  while (current != root) {
    const int next = parent[current];
    parent[current] = root;
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
    parent[root_b] = root_a;
  } else {
    parent[root_a] = root_b;
  }
}

StripeSetup BuildStripeSetup(int height, int width) {
  StripeSetup setup;

  setup.num_stripes = std::min(height, oneapi::tbb::this_task_arena::max_concurrency() * 2);
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
  tbb::parallel_for(tbb::blocked_range<int>(0, total_max_labels), [&](const tbb::blocked_range<int> &range) {
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
  if (left_label != 0 && top_label == 0) {
    labels_flat[idx] = left_label;
    return;
  }
  if (left_label == 0 && top_label != 0) {
    labels_flat[idx] = top_label;
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
                 std::vector<int> &stripe_max_used, int width, int stripe) {
  const int start_row = stripe_bounds[static_cast<std::size_t>(stripe)];
  const int end_row = stripe_bounds[static_cast<std::size_t>(stripe) + 1ULL];
  int next_label = stripe_base_label[static_cast<std::size_t>(stripe)];

  for (int row = start_row; row < end_row; ++row) {
    const auto row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width);
    const auto prev_row_offset = static_cast<std::size_t>(row - 1) * static_cast<std::size_t>(width);

    for (int col = 0; col < width; ++col) {
      const auto idx = row_offset + static_cast<std::size_t>(col);
      if (binary_flat[idx] == 0U) {
        continue;
      }

      const int left_label = (col > 0) ? labels_flat[idx - 1ULL] : 0;
      const int top_label = (row > start_row) ? labels_flat[prev_row_offset + static_cast<std::size_t>(col)] : 0;
      AssignPixelLabel(labels_flat, parent, idx, left_label, top_label, next_label);
    }
  }

  stripe_max_used[static_cast<std::size_t>(stripe)] = next_label;
}

void MergeStripeBorders(const std::vector<std::uint8_t> &binary_flat, const std::vector<int> &stripe_bounds,
                        std::vector<int> &labels_flat, std::vector<int> &parent, int width, int num_stripes) {
  for (int stripe = 0; stripe < num_stripes - 1; ++stripe) {
    const int boundary_row = stripe_bounds[static_cast<std::size_t>(stripe) + 1ULL];
    const auto row_offset = static_cast<std::size_t>(boundary_row) * static_cast<std::size_t>(width);
    const auto prev_row_offset = static_cast<std::size_t>(boundary_row - 1) * static_cast<std::size_t>(width);

    for (int col = 0; col < width; ++col) {
      const auto bottom_idx = row_offset + static_cast<std::size_t>(col);
      const auto top_idx = prev_row_offset + static_cast<std::size_t>(col);
      if (binary_flat[bottom_idx] == 1U && binary_flat[top_idx] == 1U) {
        UnionLabels(parent, labels_flat[bottom_idx], labels_flat[top_idx]);
      }
    }
  }
}

std::vector<int> BuildCompactedLabels(std::vector<int> &parent, const std::vector<int> &stripe_base_label,
                                      const std::vector<int> &stripe_max_used, int total_max_labels, int num_stripes) {
  std::vector<int> compacted(static_cast<std::size_t>(total_max_labels), 0);
  int next_compact_id = 1;

  for (int stripe = 0; stripe < num_stripes; ++stripe) {
    for (int label = stripe_base_label[static_cast<std::size_t>(stripe)];
         label < stripe_max_used[static_cast<std::size_t>(stripe)]; ++label) {
      const int root = FindRoot(parent, label);
      if (compacted[static_cast<std::size_t>(root)] == 0) {
        compacted[static_cast<std::size_t>(root)] = next_compact_id++;
      }
      compacted[static_cast<std::size_t>(label)] = compacted[static_cast<std::size_t>(root)];
    }
  }

  return compacted;
}

void ApplyCompactedLabels(std::vector<int> &labels_flat, const std::vector<int> &compacted,
                          const std::vector<int> &stripe_bounds, int width, int num_stripes) {
  tbb::parallel_for(tbb::blocked_range<int>(0, num_stripes), [&](const tbb::blocked_range<int> &range) {
    for (int stripe = range.begin(); stripe != range.end(); ++stripe) {
      const int start_row = stripe_bounds[static_cast<std::size_t>(stripe)];
      const int end_row = stripe_bounds[static_cast<std::size_t>(stripe) + 1ULL];

      for (int row = start_row; row < end_row; ++row) {
        const auto row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width);
        for (int col = 0; col < width; ++col) {
          const auto idx = row_offset + static_cast<std::size_t>(col);
          const int label = labels_flat[idx];
          if (label != 0) {
            labels_flat[idx] = compacted[static_cast<std::size_t>(label)];
          }
        }
      }
    }
  });
}

}  // namespace

MarinLMarkComponentsTBB::MarinLMarkComponentsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool MarinLMarkComponentsTBB::IsBinary(const Image &img) {
  for (const auto &row : img) {
    for (int pixel : row) {
      if (pixel != 0 && pixel != 1) {
        return false;
      }
    }
  }
  return true;
}

bool MarinLMarkComponentsTBB::ValidationImpl() {
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

bool MarinLMarkComponentsTBB::PreProcessingImpl() {
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

  binary_flat_.assign(static_cast<std::size_t>(total_pixels), 0);
  labels_flat_.assign(static_cast<std::size_t>(total_pixels), 0);
  labels_out_.clear();

  tbb::parallel_for(tbb::blocked_range<int>(0, height_), [&](const tbb::blocked_range<int> &range) {
    for (int row = range.begin(); row != range.end(); ++row) {
      const auto row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width_);
      for (int col = 0; col < width_; ++col) {
        binary_flat_[row_offset + static_cast<std::size_t>(col)] =
            static_cast<std::uint8_t>(img[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]);
      }
    }
  });

  return true;
}

bool MarinLMarkComponentsTBB::RunImpl() {
  if (height_ == 0 || width_ == 0) {
    return true;
  }

  const StripeSetup setup = BuildStripeSetup(height_, width_);
  std::vector<int> parent(static_cast<std::size_t>(setup.total_max_labels), 0);
  InitializeParents(parent, setup.total_max_labels);

  std::vector<int> stripe_max_used(static_cast<std::size_t>(setup.num_stripes), 0);
  tbb::parallel_for(tbb::blocked_range<int>(0, setup.num_stripes), [&](const tbb::blocked_range<int> &range) {
    for (int stripe = range.begin(); stripe != range.end(); ++stripe) {
      LabelStripe(binary_flat_, labels_flat_, parent, setup.stripe_bounds, setup.stripe_base_label, stripe_max_used,
                  width_, stripe);
    }
  });

  MergeStripeBorders(binary_flat_, setup.stripe_bounds, labels_flat_, parent, width_, setup.num_stripes);
  const std::vector<int> compacted =
      BuildCompactedLabels(parent, setup.stripe_base_label, stripe_max_used, setup.total_max_labels, setup.num_stripes);
  ApplyCompactedLabels(labels_flat_, compacted, setup.stripe_bounds, width_, setup.num_stripes);

  return true;
}

bool MarinLMarkComponentsTBB::PostProcessingImpl() {
  labels_out_.clear();
  labels_out_.resize(static_cast<std::size_t>(height_));

  tbb::parallel_for(tbb::blocked_range<int>(0, height_), [&](const tbb::blocked_range<int> &range) {
    for (int row = range.begin(); row != range.end(); ++row) {
      labels_out_[static_cast<std::size_t>(row)].resize(static_cast<std::size_t>(width_));
      const auto row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width_);
      std::copy(labels_flat_.begin() + static_cast<std::ptrdiff_t>(row_offset),
                labels_flat_.begin() + static_cast<std::ptrdiff_t>(row_offset) + width_,
                labels_out_[static_cast<std::size_t>(row)].begin());
    }
  });

  GetOutput().labels = std::move(labels_out_);
  return true;
}

}  // namespace marin_l_mark_components
