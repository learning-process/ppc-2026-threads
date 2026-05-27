#include "klimovich_v_crs_complex_mat_mul/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "klimovich_v_crs_complex_mat_mul/common/include/common.hpp"
#include "util/include/util.hpp"

namespace klimovich_v_crs_complex_mat_mul {
namespace {

struct RowStage {
  std::vector<int> cols;
  std::vector<Cplx> vals;
};

void GustavsonRow(const CrsMatrix &lhs, const CrsMatrix &rhs, int row, std::vector<Cplx> &spa,
                  std::vector<int> &touched_by_row, std::vector<int> &touched_cols, RowStage &stage) {
  touched_cols.clear();

  for (int lp = lhs.row_offsets[row]; lp < lhs.row_offsets[row + 1]; ++lp) {
    const int k = lhs.col_indices[lp];
    const Cplx a_ik = lhs.data[lp];
    for (int rq = rhs.row_offsets[k]; rq < rhs.row_offsets[k + 1]; ++rq) {
      const int j = rhs.col_indices[rq];
      if (touched_by_row[j] != row) {
        touched_by_row[j] = row;
        touched_cols.push_back(j);
        spa[j] = a_ik * rhs.data[rq];
      } else {
        spa[j] += a_ik * rhs.data[rq];
      }
    }
  }

  std::ranges::sort(touched_cols);

  stage.cols.clear();
  stage.vals.clear();
  stage.cols.reserve(touched_cols.size());
  stage.vals.reserve(touched_cols.size());

  for (const int j : touched_cols) {
    const Cplx v = spa[j];
    spa[j] = Cplx(0.0, 0.0);
    if (std::abs(v.real()) > kZeroDropTol || std::abs(v.imag()) > kZeroDropTol) {
      stage.cols.push_back(j);
      stage.vals.push_back(v);
    }
  }
}

CrsMatrix Assemble(int rows, int cols, const std::vector<RowStage> &per_row) {
  CrsMatrix out(rows, cols);
  for (int i = 0; i < rows; ++i) {
    out.row_offsets[i + 1] = out.row_offsets[i] + static_cast<int>(per_row[i].cols.size());
  }
  out.col_indices.reserve(static_cast<std::size_t>(out.row_offsets[rows]));
  out.data.reserve(static_cast<std::size_t>(out.row_offsets[rows]));
  for (int i = 0; i < rows; ++i) {
    out.col_indices.insert(out.col_indices.end(), per_row[i].cols.begin(), per_row[i].cols.end());
    out.data.insert(out.data.end(), per_row[i].vals.begin(), per_row[i].vals.end());
  }
  return out;
}

}  // namespace

KlimovichVCrsComplexMatMulStl::KlimovichVCrsComplexMatMulStl(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CrsMatrix();
}

bool KlimovichVCrsComplexMatMulStl::ValidationImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  return lhs.n_cols == rhs.n_rows;
}

bool KlimovichVCrsComplexMatMulStl::PreProcessingImpl() {
  return true;
}

CrsMatrix KlimovichVCrsComplexMatMulStl::MultiplyCrs(const CrsMatrix &lhs, const CrsMatrix &rhs) {
  std::vector<RowStage> per_row(static_cast<std::size_t>(lhs.n_rows));

  const int requested_threads = std::max(1, ppc::util::GetNumThreads());
  const int total_rows = lhs.n_rows;
  const int worker_count = std::max(1, std::min(requested_threads, total_rows));

  auto worker = [&](int begin, int end) {
    std::vector<Cplx> spa(static_cast<std::size_t>(rhs.n_cols), Cplx(0.0, 0.0));
    std::vector<int> touched_by_row(static_cast<std::size_t>(rhs.n_cols), -1);
    std::vector<int> touched_cols;
    touched_cols.reserve(static_cast<std::size_t>(rhs.n_cols));

    for (int i = begin; i < end; ++i) {
      GustavsonRow(lhs, rhs, i, spa, touched_by_row, touched_cols, per_row[i]);
    }
  };

  const int chunk = total_rows / worker_count;
  const int leftover = total_rows % worker_count;

  std::vector<std::thread> pool;
  pool.reserve(static_cast<std::size_t>(worker_count));

  int cursor = 0;
  for (int worker_index = 0; worker_index < worker_count; ++worker_index) {
    const int begin = cursor;
    const int end = begin + chunk + (worker_index < leftover ? 1 : 0);
    cursor = end;
    if (begin < end) {
      pool.emplace_back(worker, begin, end);
    }
  }
  for (auto &th : pool) {
    th.join();
  }

  return Assemble(lhs.n_rows, rhs.n_cols, per_row);
}

bool KlimovichVCrsComplexMatMulStl::RunImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  GetOutput() = MultiplyCrs(lhs, rhs);
  return true;
}

bool KlimovichVCrsComplexMatMulStl::PostProcessingImpl() {
  return true;
}

}  // namespace klimovich_v_crs_complex_mat_mul
