#include "kondrashova_v_marking_components_seq/seq/include/ops_seq.hpp"

#include <queue>
#include <utility>
#include <vector>

namespace kondrashova_v_marking_components_seq {

KondrashovaVTaskSEQ::KondrashovaVTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool KondrashovaVTaskSEQ::ValidationImpl() {
  // Фреймворк может вызывать Validation до установки реальных данных
  return true;
}

bool KondrashovaVTaskSEQ::PreProcessingImpl() {
  const auto &in = GetInput();

  width_ = in.width;
  height_ = in.height;
  image_ = in.data;

  if (width_ > 0 && height_ > 0 && 
      static_cast<int>(image_.size()) == width_ * height_) {
    labels_1d_.assign(width_ * height_, 0);
  } else {
    labels_1d_.clear();
  }

  GetOutput().count = 0;
  GetOutput().labels.clear();
  return true;
}

bool KondrashovaVTaskSEQ::RunImpl() {
  if (width_ <= 0 || height_ <= 0 || image_.empty()) {
    GetOutput().count = 0;
    return true;
  }

  int current_label = 0;

  auto inside = [&](int x, int y) {
    return (x >= 0 && x < height_ && y >= 0 && y < width_);
  };

  const int dx[4] = {-1, 1, 0, 0};
  const int dy[4] = {0, 0, -1, 1};

  for (int i = 0; i < height_; ++i) {
    for (int j = 0; j < width_; ++j) {
      int idx = i * width_ + j;

      if (image_[idx] == 0 && labels_1d_[idx] == 0) {
        current_label++;

        std::queue<std::pair<int, int>> q;
        q.emplace(i, j);
        labels_1d_[idx] = current_label;

        while (!q.empty()) {
          auto [x, y] = q.front();
          q.pop();

          for (int k = 0; k < 4; ++k) {
            int nx = x + dx[k];
            int ny = y + dy[k];

            if (!inside(nx, ny)) continue;

            int nidx = nx * width_ + ny;

            if (image_[nidx] == 0 && labels_1d_[nidx] == 0) {
              labels_1d_[nidx] = current_label;
              q.emplace(nx, ny);
            }
          }
        }
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

  for (int i = 0; i < height_; ++i) {
    for (int j = 0; j < width_; ++j) {
      GetOutput().labels[i][j] = labels_1d_[i * width_ + j];
    }
  }

  return true;
}

}  // namespace kondrashova_v_marking_components_seq