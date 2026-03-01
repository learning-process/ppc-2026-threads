#include "kondrashova_v_marking_components_seq/seq/include/ops_seq.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

#include "kondrashova_v_marking_components_seq/common/include/common.hpp"

namespace kondrashova_v_marking_components_seq {

KondrashovaVTaskSEQ::KondrashovaVTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool KondrashovaVTaskSEQ::ValidationImpl() {
  return true;
}

bool KondrashovaVTaskSEQ::PreProcessingImpl() {
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

void BfsComponent(int start_i, int start_j, int width, int height, int label, const std::vector<uint8_t> &image,
                  std::vector<int> &labels_1d) {
  const std::array<int, 4> dx = {-1, 1, 0, 0};
  const std::array<int, 4> dy = {0, 0, -1, 1};

  auto inside = [&](int xi, int yi) { return (xi >= 0 && xi < height && yi >= 0 && yi < width); };

  std::queue<std::pair<int, int>> q;
  q.emplace(start_i, start_j);
  labels_1d[(static_cast<size_t>(start_i) * static_cast<size_t>(width)) + static_cast<size_t>(start_j)] = label;

  while (!q.empty()) {
    auto [cx, cy] = q.front();
    q.pop();

    for (int ki = 0; ki < 4; ++ki) {
      int nx = cx + dx.at(ki);
      int ny = cy + dy.at(ki);

      if (!inside(nx, ny)) {
        continue;
      }

      auto nidx = (static_cast<size_t>(nx) * static_cast<size_t>(width)) + static_cast<size_t>(ny);
      if (image[nidx] == 0 && labels_1d[nidx] == 0) {
        labels_1d[nidx] = label;
        q.emplace(nx, ny);
      }
    }
  }
}

}  // namespace

bool KondrashovaVTaskSEQ::RunImpl() {
  if (width_ <= 0 || height_ <= 0 || image_.empty()) {
    GetOutput().count = 0;
    return true;
  }

  int current_label = 0;

  for (int ii = 0; ii < height_; ++ii) {
    for (int jj = 0; jj < width_; ++jj) {
      auto idx = (static_cast<size_t>(ii) * static_cast<size_t>(width_)) + static_cast<size_t>(jj);
      if (image_[idx] == 0 && labels_1d_[idx] == 0) {
        ++current_label;
        BfsComponent(ii, jj, width_, height_, current_label, image_, labels_1d_);
      }
    }
  }

  GetOutput().count = current_label;
  return true;
}

bool KondrashovaVTaskSEQ::PostProcessingImpl() {
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

}  // namespace kondrashova_v_marking_components_seq
