#include "ermakov_a_spar_mat_mult/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

#include "ermakov_a_spar_mat_mult/common/include/common.hpp"
#include "task/include/task.hpp"
#include "util/include/util.hpp"

namespace ermakov_a_spar_mat_mult {

namespace {

struct LocalRowData {
  std::vector<int> cols;
  std::vector<std::complex<double>> vals;
};

std::vector<int> BuildRowBounds(const MatrixCRS &a, const MatrixCRS &b, int proc_count) {
  if (proc_count <= 0) {
    return {};
  }

  std::vector<int> bounds(static_cast<std::size_t>(proc_count) + 1ULL, 0);
  bounds.back() = a.rows;

  std::vector<int> row_costs(static_cast<std::size_t>(a.rows), 0);
  std::int64_t total_cost = 0;
  for (int row = 0; row < a.rows; ++row) {
    int row_cost = 0;
    for (int ak = a.row_ptr[static_cast<std::size_t>(row)]; ak < a.row_ptr[static_cast<std::size_t>(row) + 1ULL];
         ++ak) {
      const int b_row = a.col_index[static_cast<std::size_t>(ak)];
      row_cost += b.row_ptr[static_cast<std::size_t>(b_row) + 1ULL] - b.row_ptr[static_cast<std::size_t>(b_row)];
    }
    row_costs[static_cast<std::size_t>(row)] = row_cost;
    total_cost += row_cost;
  }

  if (proc_count <= 1 || total_cost == 0) {
    for (int proc = 0; proc <= proc_count; ++proc) {
      bounds[static_cast<std::size_t>(proc)] = (proc * a.rows) / proc_count;
    }
    return bounds;
  }

  int next_proc = 1;
  std::int64_t prefix_cost = 0;
  for (int row = 0; row < a.rows && next_proc < proc_count; ++row) {
    prefix_cost += row_costs[static_cast<std::size_t>(row)];
    const std::int64_t target_cost = (static_cast<std::int64_t>(next_proc) * total_cost) / proc_count;
    if (prefix_cost >= target_cost) {
      bounds[static_cast<std::size_t>(next_proc)] = row + 1;
      ++next_proc;
    }
  }

  while (next_proc < proc_count) {
    bounds[static_cast<std::size_t>(next_proc)] = a.rows;
    ++next_proc;
  }

  return bounds;
}

std::vector<int> BuildCountsFromBounds(const std::vector<int> &bounds) {
  std::vector<int> counts(bounds.size() - 1ULL, 0);
  for (std::size_t proc = 0; proc + 1 < bounds.size(); ++proc) {
    counts[proc] = bounds[proc + 1] - bounds[proc];
  }
  return counts;
}

std::vector<int> BuildDisplacements(const std::vector<int> &counts) {
  std::vector<int> displs(counts.size(), 0);
  for (std::size_t proc = 1; proc < counts.size(); ++proc) {
    displs[proc] = displs[proc - 1] + counts[proc - 1];
  }
  return displs;
}

std::vector<int> BuildNNZCounts(const MatrixCRS &matrix, const std::vector<int> &bounds) {
  std::vector<int> counts(bounds.size() - 1ULL, 0);
  for (std::size_t proc = 0; proc + 1 < bounds.size(); ++proc) {
    counts[proc] = matrix.row_ptr[static_cast<std::size_t>(bounds[proc + 1])] -
                   matrix.row_ptr[static_cast<std::size_t>(bounds[proc])];
  }
  return counts;
}

MPI_Datatype GetComplexBytesType() {
  static MPI_Datatype datatype = MPI_DATATYPE_NULL;
  static bool initialized = false;
  if (!initialized) {
    MPI_Type_contiguous(static_cast<int>(sizeof(std::complex<double>)), MPI_BYTE, &datatype);
    MPI_Type_commit(&datatype);
    initialized = true;
  }
  return datatype;
}

MatrixCRS ScatterRows(const MatrixCRS &matrix, const std::vector<int> &row_bounds, const std::vector<int> &nnz_counts,
                      int rank) {
  const std::vector<int> row_counts = BuildCountsFromBounds(row_bounds);
  const std::vector<int> row_displs = BuildDisplacements(row_counts);
  const std::vector<int> nnz_displs = BuildDisplacements(nnz_counts);

  MatrixCRS local;
  local.rows = row_counts[static_cast<std::size_t>(rank)];
  local.cols = matrix.cols;
  local.row_ptr.assign(static_cast<std::size_t>(local.rows) + 1ULL, 0);
  local.col_index.resize(static_cast<std::size_t>(nnz_counts[static_cast<std::size_t>(rank)]));
  local.values.resize(static_cast<std::size_t>(nnz_counts[static_cast<std::size_t>(rank)]));

  std::vector<int> all_row_lengths;

  if (rank == 0) {
    all_row_lengths.resize(static_cast<std::size_t>(matrix.rows), 0);
    for (int row = 0; row < matrix.rows; ++row) {
      all_row_lengths[static_cast<std::size_t>(row)] =
          matrix.row_ptr[static_cast<std::size_t>(row) + 1ULL] - matrix.row_ptr[static_cast<std::size_t>(row)];
    }
  }

  std::vector<int> local_row_lengths(static_cast<std::size_t>(local.rows), 0);
  MPI_Scatterv(all_row_lengths.data(), row_counts.data(), row_displs.data(), MPI_INT, local_row_lengths.data(),
               local.rows, MPI_INT, 0, MPI_COMM_WORLD);

  const int local_nnz = nnz_counts[static_cast<std::size_t>(rank)];
  MPI_Scatterv(matrix.col_index.data(), nnz_counts.data(), nnz_displs.data(), MPI_INT, local.col_index.data(),
               local_nnz, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatterv(matrix.values.data(), nnz_counts.data(), nnz_displs.data(), GetComplexBytesType(), local.values.data(),
               local_nnz, GetComplexBytesType(), 0, MPI_COMM_WORLD);

  int prefix = 0;
  for (int row = 0; row < local.rows; ++row) {
    local.row_ptr[static_cast<std::size_t>(row)] = prefix;
    prefix += local_row_lengths[static_cast<std::size_t>(row)];
  }
  local.row_ptr[static_cast<std::size_t>(local.rows)] = prefix;

  return local;
}

void BroadcastMatrix(MatrixCRS &matrix, int rank) {
  std::array<int, 3> dims = {matrix.rows, matrix.cols, static_cast<int>(matrix.values.size())};
  MPI_Bcast(dims.data(), static_cast<int>(dims.size()), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    matrix.rows = dims[0];
    matrix.cols = dims[1];
    matrix.values.resize(static_cast<std::size_t>(dims[2]));
    matrix.col_index.resize(static_cast<std::size_t>(dims[2]));
    matrix.row_ptr.resize(static_cast<std::size_t>(matrix.rows) + 1ULL);
  }

  MPI_Bcast(matrix.col_index.data(), dims[2], MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(matrix.row_ptr.data(), matrix.rows + 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(matrix.values.data(), dims[2], GetComplexBytesType(), 0, MPI_COMM_WORLD);
}

void AccumulateRowProducts(const MatrixCRS &a, const MatrixCRS &b, int row_index,
                           std::vector<std::complex<double>> &row_vals, std::vector<int> &row_mark,
                           std::vector<int> &used_cols) {
  used_cols.clear();

  for (int ak = a.row_ptr[static_cast<std::size_t>(row_index)];
       ak < a.row_ptr[static_cast<std::size_t>(row_index) + 1ULL]; ++ak) {
    const int b_row = a.col_index[static_cast<std::size_t>(ak)];
    const auto a_val = a.values[static_cast<std::size_t>(ak)];

    for (int bk = b.row_ptr[static_cast<std::size_t>(b_row)]; bk < b.row_ptr[static_cast<std::size_t>(b_row) + 1ULL];
         ++bk) {
      const int col = b.col_index[static_cast<std::size_t>(bk)];
      const auto product = a_val * b.values[static_cast<std::size_t>(bk)];

      if (row_mark[static_cast<std::size_t>(col)] != row_index) {
        row_mark[static_cast<std::size_t>(col)] = row_index;
        row_vals[static_cast<std::size_t>(col)] = product;
        used_cols.push_back(col);
      } else {
        row_vals[static_cast<std::size_t>(col)] += product;
      }
    }
  }
}

void CollectRowValues(const std::vector<std::complex<double>> &row_vals, std::vector<int> &used_cols,
                      LocalRowData &row) {
  std::ranges::sort(used_cols);
  row.cols.clear();
  row.vals.clear();
  row.cols.reserve(used_cols.size());
  row.vals.reserve(used_cols.size());

  for (int col : used_cols) {
    const auto &value = row_vals[static_cast<std::size_t>(col)];
    if (value != std::complex<double>(0.0, 0.0)) {
      row.cols.push_back(col);
      row.vals.push_back(value);
    }
  }
}

MatrixCRS MultiplyLocalOMP(const MatrixCRS &a, const MatrixCRS &b) {
  MatrixCRS result;
  result.rows = a.rows;
  result.cols = b.cols;
  result.row_ptr.assign(static_cast<std::size_t>(result.rows) + 1ULL, 0);

  if (a.rows == 0 || b.cols == 0) {
    return result;
  }

  int rank_count = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &rank_count);

  int thread_count = ppc::util::GetNumThreads();
  const unsigned hw_threads = std::thread::hardware_concurrency();
  if (hw_threads > 0U && rank_count > 1) {
    const unsigned per_rank_cap = std::max(1U, hw_threads / static_cast<unsigned>(rank_count));
    thread_count = std::min(thread_count, static_cast<int>(per_rank_cap));
  }
  thread_count = std::max(1, std::min(thread_count, a.rows));
  std::vector<LocalRowData> rows_data(static_cast<std::size_t>(a.rows));

#pragma omp parallel default(none) shared(a, b, rows_data) num_threads(thread_count) if (thread_count > 1)
  {
    std::vector<std::complex<double>> row_vals(static_cast<std::size_t>(b.cols), std::complex<double>(0.0, 0.0));
    std::vector<int> row_mark(static_cast<std::size_t>(b.cols), -1);
    std::vector<int> used_cols;
    used_cols.reserve(256);

#pragma omp for
    for (int row = 0; row < a.rows; ++row) {
      AccumulateRowProducts(a, b, row, row_vals, row_mark, used_cols);
      CollectRowValues(row_vals, used_cols, rows_data[static_cast<std::size_t>(row)]);
    }
  }

  int total_nnz = 0;
  for (int row = 0; row < result.rows; ++row) {
    result.row_ptr[static_cast<std::size_t>(row)] = total_nnz;
    total_nnz += static_cast<int>(rows_data[static_cast<std::size_t>(row)].vals.size());
  }
  result.row_ptr[static_cast<std::size_t>(result.rows)] = total_nnz;

  result.values.reserve(static_cast<std::size_t>(total_nnz));
  result.col_index.reserve(static_cast<std::size_t>(total_nnz));

  for (int row = 0; row < result.rows; ++row) {
    const auto &row_data = rows_data[static_cast<std::size_t>(row)];
    result.col_index.insert(result.col_index.end(), row_data.cols.begin(), row_data.cols.end());
    result.values.insert(result.values.end(), row_data.vals.begin(), row_data.vals.end());
  }

  return result;
}

void GatherMatrix(const MatrixCRS &local, MatrixCRS &global, const std::vector<int> &row_bounds, int rank, int size,
                  int total_rows) {
  const std::vector<int> row_counts = BuildCountsFromBounds(row_bounds);
  const std::vector<int> row_displs = BuildDisplacements(row_counts);
  std::vector<int> nnz_counts(static_cast<std::size_t>(size), 0);
  const int local_nnz = static_cast<int>(local.values.size());
  MPI_Gather(&local_nnz, 1, MPI_INT, nnz_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> local_row_lengths(static_cast<std::size_t>(local.rows), 0);
  for (int row = 0; row < local.rows; ++row) {
    local_row_lengths[static_cast<std::size_t>(row)] =
        local.row_ptr[static_cast<std::size_t>(row) + 1ULL] - local.row_ptr[static_cast<std::size_t>(row)];
  }

  std::vector<int> nnz_displs;
  std::vector<int> gathered_row_lengths;
  std::vector<int> gathered_cols;
  std::vector<std::complex<double>> gathered_values;

  if (rank == 0) {
    nnz_displs.resize(static_cast<std::size_t>(size), 0);
    for (int proc = 1; proc < size; ++proc) {
      nnz_displs[static_cast<std::size_t>(proc)] =
          nnz_displs[static_cast<std::size_t>(proc - 1)] + nnz_counts[static_cast<std::size_t>(proc - 1)];
    }

    gathered_row_lengths.resize(static_cast<std::size_t>(total_rows), 0);
    gathered_cols.resize(static_cast<std::size_t>(std::accumulate(nnz_counts.begin(), nnz_counts.end(), 0)), 0);
    gathered_values.resize(gathered_cols.size());
  }

  MPI_Gatherv(local_row_lengths.data(), local.rows, MPI_INT, gathered_row_lengths.data(), row_counts.data(),
              row_displs.data(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gatherv(local.col_index.data(), local_nnz, MPI_INT, gathered_cols.data(), nnz_counts.data(), nnz_displs.data(),
              MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Gatherv(local.values.data(), local_nnz, GetComplexBytesType(), gathered_values.data(), nnz_counts.data(),
              nnz_displs.data(), GetComplexBytesType(), 0, MPI_COMM_WORLD);

  if (rank != 0) {
    return;
  }

  global.row_ptr.assign(static_cast<std::size_t>(total_rows) + 1ULL, 0);
  int prefix = 0;
  for (int row = 0; row < total_rows; ++row) {
    global.row_ptr[static_cast<std::size_t>(row)] = prefix;
    prefix += gathered_row_lengths[static_cast<std::size_t>(row)];
  }
  global.row_ptr[static_cast<std::size_t>(total_rows)] = prefix;

  global.col_index = std::move(gathered_cols);
  global.values = std::move(gathered_values);
}

}  // namespace

ErmakovASparMatMultALL::ErmakovASparMatMultALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ErmakovASparMatMultALL::ValidateMatrix(const MatrixCRS &m) {
  if (m.rows < 0 || m.cols < 0) {
    return false;
  }
  if (m.row_ptr.size() != static_cast<std::size_t>(m.rows) + 1ULL) {
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

  for (int row = 0; row < m.rows; ++row) {
    if (m.row_ptr[static_cast<std::size_t>(row)] > m.row_ptr[static_cast<std::size_t>(row) + 1ULL]) {
      return false;
    }
  }

  for (int idx = 0; idx < nnz; ++idx) {
    if (m.col_index[static_cast<std::size_t>(idx)] < 0 || m.col_index[static_cast<std::size_t>(idx)] >= m.cols) {
      return false;
    }
  }

  return true;
}

bool ErmakovASparMatMultALL::ValidationImpl() {
  const auto &a = GetInput().A;
  const auto &b = GetInput().B;
  return a.cols == b.rows && ValidateMatrix(a) && ValidateMatrix(b);
}

bool ErmakovASparMatMultALL::PreProcessingImpl() {
  a_ = GetInput().A;
  b_ = GetInput().B;
  c_.rows = a_.rows;
  c_.cols = b_.cols;
  c_.values.clear();
  c_.col_index.clear();
  c_.row_ptr.assign(static_cast<std::size_t>(c_.rows) + 1ULL, 0);
  return true;
}

bool ErmakovASparMatMultALL::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0 && a_.cols != b_.rows) {
    return false;
  }

  if (size == 1) {
    c_ = MultiplyLocalOMP(a_, b_);
    return true;
  }

  BroadcastMatrix(b_, rank);

  std::vector<int> row_bounds(static_cast<std::size_t>(size) + 1ULL, 0);
  std::vector<int> nnz_counts(static_cast<std::size_t>(size), 0);
  if (rank == 0) {
    row_bounds = BuildRowBounds(a_, b_, size);
    nnz_counts = BuildNNZCounts(a_, row_bounds);
  }

  MPI_Bcast(row_bounds.data(), size + 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nnz_counts.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

  const MatrixCRS local_a = ScatterRows(a_, row_bounds, nnz_counts, rank);
  const MatrixCRS local_c = MultiplyLocalOMP(local_a, b_);

  c_.rows = a_.rows;
  c_.cols = b_.cols;
  c_.values.clear();
  c_.col_index.clear();
  c_.row_ptr.assign(static_cast<std::size_t>(c_.rows) + 1ULL, 0);
  GatherMatrix(local_c, c_, row_bounds, rank, size, a_.rows);

  if (GetStateOfTesting() == ppc::task::StateOfTesting::kPerf) {
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
  }

  BroadcastMatrix(c_, rank);
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool ErmakovASparMatMultALL::PostProcessingImpl() {
  GetOutput() = c_;
  return true;
}

}  // namespace ermakov_a_spar_mat_mult
