#include "marin_l_mark_components/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#include "marin_l_mark_components/common/include/common.hpp"
#include "util/include/util.hpp"

namespace marin_l_mark_components {

namespace {

constexpr std::uint64_t kMaxPixels = 100000000ULL;
constexpr int kMinRowsPerStripe = 64;

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

int MergeLabels(std::vector<int> &parent, int left_label, int top_label, int &next_label) {
  if (left_label == 0 && top_label == 0) {
    parent[static_cast<std::size_t>(next_label)] = next_label;
    return next_label++;
  }

  if (left_label == 0) {
    return top_label;
  }

  if (top_label == 0) {
    return left_label;
  }

  const int merged_label = std::min(left_label, top_label);
  if (left_label != top_label) {
    UnionLabels(parent, left_label, top_label);
  }

  return merged_label;
}

void RunParallelWorkers(int worker_count, const std::function<void(int, int)> &fn) {
  const int total_workers = std::max(1, worker_count);
  std::vector<std::thread> threads;
  threads.reserve(static_cast<std::size_t>(std::max(0, total_workers - 1)));

  std::exception_ptr thread_error;
  std::mutex error_mutex;

  for (int worker = 1; worker < total_workers; ++worker) {
    threads.emplace_back([&, worker]() {
      try {
        fn(worker, total_workers);
      } catch (...) {
        std::scoped_lock lock(error_mutex);
        if (thread_error == nullptr) {
          thread_error = std::current_exception();
        }
      }
    });
  }

  try {
    fn(0, total_workers);
  } catch (...) {
    std::scoped_lock lock(error_mutex);
    if (thread_error == nullptr) {
      thread_error = std::current_exception();
    }
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  if (thread_error != nullptr) {
    std::rethrow_exception(thread_error);
  }
}

void ParallelForBlocks(int total_items, int worker_count, const std::function<void(int, int)> &fn) {
  if (total_items <= 0) {
    return;
  }

  const int actual_workers = std::min(std::max(1, worker_count), total_items);
  RunParallelWorkers(actual_workers, [fn, total_items](int worker, int total_workers) {
    const int begin = (worker * total_items) / total_workers;
    const int end = ((worker + 1) * total_items) / total_workers;
    if (begin < end) {
      fn(begin, end);
    }
  });
}

void InitializeParents(std::vector<int> &parent, int total_max_labels, int worker_count) {
  ParallelForBlocks(total_max_labels, worker_count, [&](int begin, int end) {
    for (int label = begin; label < end; ++label) {
      parent[static_cast<std::size_t>(label)] = label;
    }
  });
}

void ProcessStripe(const std::vector<std::uint8_t> &binary, std::vector<int> &labels_flat, std::vector<int> &parent,
                   const std::vector<int> &stripe_bounds, const std::vector<int> &stripe_base_label,
                   std::vector<int> &stripe_used_label_end, int width, int stripe) {
  const int row_begin = stripe_bounds[static_cast<std::size_t>(stripe)];
  const int row_end = stripe_bounds[static_cast<std::size_t>(stripe) + 1ULL];
  int next_label = stripe_base_label[static_cast<std::size_t>(stripe)];

  for (int row = row_begin; row < row_end; ++row) {
    const std::size_t row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width);
    const std::size_t top_row_offset = row_offset - static_cast<std::size_t>(width);

    for (int col = 0; col < width; ++col) {
      const std::size_t idx = row_offset + static_cast<std::size_t>(col);
      if (binary[idx] == 0U) {
        continue;
      }

      const int left_label = (col > 0) ? labels_flat[idx - 1ULL] : 0;
      const int top_label = (row > row_begin) ? labels_flat[top_row_offset + static_cast<std::size_t>(col)] : 0;
      labels_flat[idx] = MergeLabels(parent, left_label, top_label, next_label);
    }
  }

  stripe_used_label_end[static_cast<std::size_t>(stripe)] = next_label;
}

std::vector<int> BuildCompactedLabels(std::vector<int> &parent, const std::vector<int> &stripe_base_label,
                                      const std::vector<int> &stripe_used_label_end) {
  std::vector<int> compacted(parent.size(), 0);
  int next_compact_id = 1;

  for (std::size_t stripe = 0; stripe < stripe_base_label.size(); ++stripe) {
    for (int label = stripe_base_label[stripe]; label < stripe_used_label_end[stripe]; ++label) {
      const int root = FindRoot(parent, label);
      if (compacted[static_cast<std::size_t>(root)] == 0) {
        compacted[static_cast<std::size_t>(root)] = next_compact_id++;
      }
      compacted[static_cast<std::size_t>(label)] = compacted[static_cast<std::size_t>(root)];
    }
  }

  return compacted;
}

}  // namespace

MarinLMarkComponentsSTL::MarinLMarkComponentsSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool MarinLMarkComponentsSTL::IsBinary(const Image &img) {
  for (const auto &row : img) {
    for (int pixel : row) {
      if (pixel != 0 && pixel != 1) {
        return false;
      }
    }
  }
  return true;
}

bool MarinLMarkComponentsSTL::ValidationImpl() {
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

bool MarinLMarkComponentsSTL::PreProcessingImpl() {
  const auto &input_binary = GetInput().binary;
  height_ = static_cast<int>(input_binary.size());
  width_ = static_cast<int>(input_binary.front().size());

  if (height_ <= 0 || width_ <= 0) {
    return false;
  }

  const std::uint64_t pixels = static_cast<std::uint64_t>(height_) * static_cast<std::uint64_t>(width_);
  if (pixels > kMaxPixels) {
    return false;
  }

  const int requested_threads = std::max(1, ppc::util::GetNumThreads());
  stripe_count_ = std::min(height_, requested_threads * 2);
  if (stripe_count_ > 0 && height_ / stripe_count_ < kMinRowsPerStripe) {
    stripe_count_ = std::max(1, height_ / kMinRowsPerStripe);
  }
  stripe_count_ = std::max(1, stripe_count_);

  stripe_bounds_.assign(static_cast<std::size_t>(stripe_count_) + 1ULL, 0);
  for (int stripe = 0; stripe <= stripe_count_; ++stripe) {
    stripe_bounds_[static_cast<std::size_t>(stripe)] = (stripe * height_) / stripe_count_;
  }

  stripe_base_label_.assign(static_cast<std::size_t>(stripe_count_), 0);
  stripe_used_label_end_.assign(static_cast<std::size_t>(stripe_count_), 0);
  total_max_labels_ = 1;

  for (int stripe = 0; stripe < stripe_count_; ++stripe) {
    stripe_base_label_[static_cast<std::size_t>(stripe)] = total_max_labels_;
    const int stripe_height =
        stripe_bounds_[static_cast<std::size_t>(stripe) + 1ULL] - stripe_bounds_[static_cast<std::size_t>(stripe)];
    total_max_labels_ += ((stripe_height * width_) / 2) + 1;
  }

  binary_.assign(static_cast<std::size_t>(pixels), 0);
  labels_flat_.assign(static_cast<std::size_t>(pixels), 0);
  labels_.clear();
  parent_.assign(static_cast<std::size_t>(total_max_labels_), 0);
  root_to_compact_.assign(static_cast<std::size_t>(total_max_labels_), 0);

  InitializeParents(parent_, total_max_labels_, requested_threads);

  ParallelForBlocks(height_, requested_threads, [&](int row_begin, int row_end) {
    for (int row = row_begin; row < row_end; ++row) {
      const std::size_t row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width_);
      for (int col = 0; col < width_; ++col) {
        binary_[row_offset + static_cast<std::size_t>(col)] =
            static_cast<std::uint8_t>(input_binary[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]);
      }
    }
  });

  return true;
}

bool MarinLMarkComponentsSTL::RunImpl() {
  if (height_ == 0 || width_ == 0) {
    return true;
  }

  FirstPassSTL();
  MergeStripeBorders();
  SecondPassSTL();
  return true;
}

void MarinLMarkComponentsSTL::FirstPassSTL() {
  std::ranges::fill(labels_flat_, 0);

  ParallelForBlocks(stripe_count_, stripe_count_, [&](int stripe_begin, int stripe_end) {
    for (int stripe = stripe_begin; stripe < stripe_end; ++stripe) {
      ProcessStripe(binary_, labels_flat_, parent_, stripe_bounds_, stripe_base_label_, stripe_used_label_end_, width_,
                    stripe);
    }
  });
}

void MarinLMarkComponentsSTL::MergeStripeBorders() {
  for (int stripe = 0; stripe < stripe_count_ - 1; ++stripe) {
    const int boundary_row = stripe_bounds_[static_cast<std::size_t>(stripe) + 1ULL];
    const std::size_t top_row_offset = static_cast<std::size_t>(boundary_row - 1) * static_cast<std::size_t>(width_);
    const std::size_t bottom_row_offset = static_cast<std::size_t>(boundary_row) * static_cast<std::size_t>(width_);

    for (int col = 0; col < width_; ++col) {
      const std::size_t top_idx = top_row_offset + static_cast<std::size_t>(col);
      const std::size_t bottom_idx = bottom_row_offset + static_cast<std::size_t>(col);

      if (binary_[top_idx] == 1U && binary_[bottom_idx] == 1U) {
        UnionLabels(parent_, labels_flat_[top_idx], labels_flat_[bottom_idx]);
      }
    }
  }
}

void MarinLMarkComponentsSTL::SecondPassSTL() {
  if (total_max_labels_ <= 1) {
    return;
  }

  root_to_compact_ = BuildCompactedLabels(parent_, stripe_base_label_, stripe_used_label_end_);

  ParallelForBlocks(static_cast<int>(labels_flat_.size()), stripe_count_, [&](int pixel_begin, int pixel_end) {
    for (int idx = pixel_begin; idx < pixel_end; ++idx) {
      const int label = labels_flat_[static_cast<std::size_t>(idx)];
      if (label == 0) {
        continue;
      }

      labels_flat_[static_cast<std::size_t>(idx)] = root_to_compact_[static_cast<std::size_t>(label)];
    }
  });
}

bool MarinLMarkComponentsSTL::PostProcessingImpl() {
  ConvertLabelsToOutput();
  OutType out;
  out.labels = labels_;
  GetOutput() = out;
  return true;
}

void MarinLMarkComponentsSTL::ConvertLabelsToOutput() {
  labels_.clear();
  labels_.resize(static_cast<std::size_t>(height_));

  ParallelForBlocks(height_, stripe_count_, [&](int row_begin, int row_end) {
    for (int row = row_begin; row < row_end; ++row) {
      labels_[static_cast<std::size_t>(row)].resize(static_cast<std::size_t>(width_));
      const std::size_t row_offset = static_cast<std::size_t>(row) * static_cast<std::size_t>(width_);
      for (int col = 0; col < width_; ++col) {
        labels_[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] =
            labels_flat_[row_offset + static_cast<std::size_t>(col)];
      }
    }
  });
}

}  // namespace marin_l_mark_components
