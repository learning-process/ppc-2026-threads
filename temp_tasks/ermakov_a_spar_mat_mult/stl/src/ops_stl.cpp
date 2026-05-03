#include "ermakov_a_spar_mat_mult/stl/include/ops_stl.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <numeric>
#include <thread>
#include <vector>

#include "ermakov_a_spar_mat_mult/common/include/common.hpp"
#include "util/include/util.hpp"

namespace ermakov_a_spar_mat_mult {

namespace {

constexpr std::complex<double> kZero{0.0, 0.0};
constexpr std::size_t kMinWorkPerThread = 4096;

}  // namespace

ErmakovASparMatMultSTL::ErmakovASparMatMultSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ErmakovASparMatMultSTL::ValidateMatrix(const MatrixCRS &m) {
  if (m.rows < 0 || m.cols < 0) {
    return false;
  }

  if (m.row_ptr.size() != static_cast<std::size_t>(m.rows) + 1) {
    return false;
  }

  if (m.values.size() != m.col_index.size()) {
    return false;
  }

  if (m.row_ptr.empty()) {
    return false;
  }

  const int nnz = static_cast<int>(m.values.size());
  if (m.row_ptr.front() != 0 || m.row_ptr.back() != nnz) {
    return false;
  }

  for (int i = 0; i < m.rows; ++i) {
    if (m.row_ptr[i] > m.row_ptr[i + 1]) {
      return false;
    }
  }

  for (int k = 0; k < nnz; ++k) {
    if (m.col_index[k] < 0 || m.col_index[k] >= m.cols) {
      return false;
    }
  }

  return true;
}

bool ErmakovASparMatMultSTL::ValidationImpl() {
  const auto &a = GetInput().A;
  const auto &b = GetInput().B;

  return a.cols == b.rows && ValidateMatrix(a) && ValidateMatrix(b);
}

bool ErmakovASparMatMultSTL::PreProcessingImpl() {
  a_ = GetInput().A;
  b_ = GetInput().B;

  c_.rows = a_.rows;
  c_.cols = b_.cols;
  c_.values.clear();
  c_.col_index.clear();
  c_.row_ptr.assign(static_cast<std::size_t>(c_.rows) + 1, 0);

  return true;
}

int ErmakovASparMatMultSTL::ResolveThreadCount(int row_count, std::size_t total_work) {
  if (row_count <= 1 || total_work < kMinWorkPerThread) {
    return 1;
  }

  const int hw_threads = std::max(1, ppc::util::GetNumThreads());
  const std::size_t work_limited_threads = std::max<std::size_t>(1, total_work / kMinWorkPerThread);
  return std::max(1, std::min({hw_threads, row_count, static_cast<int>(work_limited_threads)}));
}

std::vector<int> ErmakovASparMatMultSTL::BuildRowCosts(const MatrixCRS &a, const MatrixCRS &b) {
  std::vector<int> row_costs(static_cast<std::size_t>(a.rows), 1);

  for (int i = 0; i < a.rows; ++i) {
    int cost = 0;
    for (int ak = a.row_ptr[i]; ak < a.row_ptr[i + 1]; ++ak) {
      const int b_row = a.col_index[ak];
      cost += b.row_ptr[b_row + 1] - b.row_ptr[b_row];
    }
    row_costs[static_cast<std::size_t>(i)] = std::max(1, cost);
  }

  return row_costs;
}

std::vector<int> ErmakovASparMatMultSTL::BuildThreadBoundaries(const std::vector<int> &row_costs, int thread_count) {
  const int safe_thread_count = std::max(1, thread_count);
  std::vector<int> boundaries(static_cast<std::size_t>(safe_thread_count) + 1, 0);
  boundaries[static_cast<std::size_t>(safe_thread_count)] = static_cast<int>(row_costs.size());

  if (safe_thread_count <= 1 || row_costs.empty()) {
    return boundaries;
  }

  const std::size_t total_work = std::accumulate(row_costs.begin(), row_costs.end(), std::size_t{0});
  const std::size_t target_chunk = std::max<std::size_t>(1, total_work / static_cast<std::size_t>(safe_thread_count));

  std::size_t accumulated = 0;
  std::size_t next_target = target_chunk;
  int boundary_index = 1;
  const auto row_count = row_costs.size();

  for (std::size_t row = 0; row < row_count && boundary_index < safe_thread_count; ++row) {
    accumulated += static_cast<std::size_t>(row_costs[row]);
    if (accumulated >= next_target) {
      boundaries[static_cast<std::size_t>(boundary_index)] = static_cast<int>(row) + 1;
      ++boundary_index;
      next_target = target_chunk * static_cast<std::size_t>(boundary_index);
    }
  }

  for (; boundary_index < safe_thread_count; ++boundary_index) {
    boundaries[static_cast<std::size_t>(boundary_index)] = boundaries[static_cast<std::size_t>(boundary_index) - 1];
  }

  return boundaries;
}

void ErmakovASparMatMultSTL::ResetWorkspace(Workspace &workspace, int cols) {
  if (workspace.accum.size() != static_cast<std::size_t>(cols)) {
    workspace.accum.assign(static_cast<std::size_t>(cols), kZero);
    workspace.marks.assign(static_cast<std::size_t>(cols), -1);
    workspace.touched_cols.clear();
    workspace.touched_cols.reserve(256);
  }
}

void ErmakovASparMatMultSTL::MultiplyRow(int row_index, Workspace &workspace, RowData &row_data) const {
  workspace.touched_cols.clear();

  for (int ak = a_.row_ptr[row_index]; ak < a_.row_ptr[row_index + 1]; ++ak) {
    const int b_row = a_.col_index[ak];
    const auto a_value = a_.values[ak];

    for (int bk = b_.row_ptr[b_row]; bk < b_.row_ptr[b_row + 1]; ++bk) {
      const int col = b_.col_index[bk];
      const auto product = a_value * b_.values[bk];

      if (workspace.marks[static_cast<std::size_t>(col)] != row_index) {
        workspace.marks[static_cast<std::size_t>(col)] = row_index;
        workspace.accum[static_cast<std::size_t>(col)] = product;
        workspace.touched_cols.push_back(col);
      } else {
        workspace.accum[static_cast<std::size_t>(col)] += product;
      }
    }
  }

  std::ranges::sort(workspace.touched_cols);

  row_data.cols.clear();
  row_data.vals.clear();
  row_data.cols.reserve(workspace.touched_cols.size());
  row_data.vals.reserve(workspace.touched_cols.size());

  for (int col : workspace.touched_cols) {
    const auto &value = workspace.accum[static_cast<std::size_t>(col)];
    if (value != kZero) {
      row_data.cols.push_back(col);
      row_data.vals.push_back(value);
    }
  }
}

void ErmakovASparMatMultSTL::FinalizeResult(const std::vector<RowData> &rows_data) {
  int total_nnz = 0;

  for (int i = 0; i < c_.rows; ++i) {
    c_.row_ptr[static_cast<std::size_t>(i)] = total_nnz;
    total_nnz += static_cast<int>(rows_data[static_cast<std::size_t>(i)].vals.size());
  }

  c_.row_ptr[static_cast<std::size_t>(c_.rows)] = total_nnz;
  c_.values.resize(static_cast<std::size_t>(total_nnz));
  c_.col_index.resize(static_cast<std::size_t>(total_nnz));

  for (int i = 0; i < c_.rows; ++i) {
    const auto &row = rows_data[static_cast<std::size_t>(i)];
    const auto offset = static_cast<std::size_t>(c_.row_ptr[static_cast<std::size_t>(i)]);

    std::ranges::copy(row.cols, c_.col_index.begin() + static_cast<std::ptrdiff_t>(offset));
    std::ranges::copy(row.vals, c_.values.begin() + static_cast<std::ptrdiff_t>(offset));
  }
}

bool ErmakovASparMatMultSTL::RunImpl() {
  if (a_.cols != b_.rows) {
    return false;
  }

  c_.values.clear();
  c_.col_index.clear();
  std::ranges::fill(c_.row_ptr, 0);

  if (a_.rows == 0 || b_.cols == 0) {
    return true;
  }

  const std::vector<int> row_costs = BuildRowCosts(a_, b_);
  const std::size_t total_work = std::accumulate(row_costs.begin(), row_costs.end(), std::size_t{0});
  const int thread_count = ResolveThreadCount(a_.rows, total_work);

  std::vector<RowData> rows_data(static_cast<std::size_t>(a_.rows));

  if (thread_count == 1) {
    Workspace workspace;
    ResetWorkspace(workspace, b_.cols);
    for (int row = 0; row < a_.rows; ++row) {
      MultiplyRow(row, workspace, rows_data[static_cast<std::size_t>(row)]);
    }
    FinalizeResult(rows_data);
    return true;
  }

  const std::vector<int> boundaries = BuildThreadBoundaries(row_costs, thread_count);
  std::vector<std::thread> workers;
  workers.reserve(static_cast<std::size_t>(thread_count));

  for (int thread_id = 0; thread_id < thread_count; ++thread_id) {
    workers.emplace_back([&, thread_id] {
      Workspace workspace;
      ResetWorkspace(workspace, b_.cols);

      const int row_begin = boundaries[static_cast<std::size_t>(thread_id)];
      const int row_end = boundaries[static_cast<std::size_t>(thread_id) + 1];

      for (int row = row_begin; row < row_end; ++row) {
        MultiplyRow(row, workspace, rows_data[static_cast<std::size_t>(row)]);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }

  FinalizeResult(rows_data);
  return true;
}

bool ErmakovASparMatMultSTL::PostProcessingImpl() {
  GetOutput() = c_;
  return true;
}

}  // namespace ermakov_a_spar_mat_mult
